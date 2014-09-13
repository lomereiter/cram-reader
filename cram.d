import std.stdio, std.stream, std.bitmanip, std.traits, std.exception, std.conv,
       std.algorithm, std.array, std.string, std.zlib, std.range, std.typecons,
       std.parallelism, std.mmfile;

import bio.bam.abstractreader, bio.bam.read, bio.bam.referenceinfo,
    bio.bam.reference, bio.sam.header, bio.bam.tagvalue, bio.core.utils.stream,
    bio.core.base;

bool seekCur(ref const(ubyte)[] stream, size_t n) {
    if (stream.length < n)
	return false;
    stream = stream[n .. $];
    return true;
}

bool seekCur(std.stream.Stream stream, size_t n) {
    try {
	stream.seekCur(n);
    } catch(std.stream.SeekException e) {
	return false;
    }
    return true;
}

bool read(ref const(ubyte)[] stream, ubyte[] data) {
    if (stream.length < data.length)
	return false;
    auto n = data.length;
    data[0 .. n] = stream[0 .. n];
    stream = stream[n .. $];
    return true;
}

bool read(std.stream.Stream stream, ubyte[] data) {
    int _read = 0;
    auto ptr = data.ptr;
    while (_read < data.length) {
	auto n = stream.readBlock(ptr, data.length - _read);
	if (n == 0) return false;
	_read += n;
	ptr += n;
    }
    return _read == data.length;
}

bool read(S, T)(auto ref S stream, T* k)
    if (isIntegral!T)
{
    // little-endian format is assumed
    ubyte* ptr = cast(ubyte*)(k);
    ubyte[] integer_data = ptr[0 .. T.sizeof];
    return read(stream, integer_data);
}

struct itf8 {
    int value;
    alias value this;
}

bool read(S)(auto ref S stream, itf8 *k) {
    ubyte first_byte = void;
    if (!read(stream, &first_byte))
	return false;
    ubyte mask = 128;
    ubyte n;
    while (first_byte & mask) {
	++n;
	mask >>= 1;
	if (n == 4)
	    break;
    }
    k.value = first_byte & (mask - 1);
    while (n > 0) {
	if (!read(stream, &mask))
	    return false;
	--n;
	k.value <<= 8;
	k.value |= mask; // FIXME: signed?
    }
    return true;
}

bool read(S)(auto ref S stream, out itf8[] buffer) {
    itf8 length;
    if (!read(stream, &length))
	return false;
    buffer.length = length;
    foreach (i; 0 .. length)
	if (!read(stream, &buffer[i]))
	    return false;
    return true;
}

class CramReadError : Exception {
    this(string msg) {
	super(msg);
    }
}

struct CramFileDefinition {
    ubyte[4] magic;
    ubyte major;
    ubyte minor;
    ubyte[20] id;
}

bool read(S)(auto ref S stream, out CramFileDefinition def) {
    return read(stream, def.magic[]) &&
	read(stream, &def.major) &&
	read(stream, &def.minor) &&
	read(stream, def.id[]);
}

struct CramBlock {
    static enum CompressionMethod : ubyte {
	raw = 0,
	gzip = 1,
	bzip2 = 2
    }

    static enum ContentType : ubyte {
	fileHeader = 0,
	compressionHeader = 1,
	mappedSliceHeader = 2,
	externalData = 4,
	coreData = 5
    }

    CompressionMethod compression_method;
    ContentType content_type;
    itf8 content_id;
    itf8 compressed_size;
    itf8 uncompressed_size;
    const(ubyte)[] compressed_data;
    private ubyte[] uncompressed_data_;

    const(ubyte)[] uncompressed_data() {
	if (compression_method == CompressionMethod.raw)
	    return compressed_data;
	if (compression_method == CompressionMethod.gzip) {
	    if (uncompressed_data_.length != uncompressed_size) {
		auto data = cast(ubyte[])(compressed_data);
		uncompressed_data_ = cast(ubyte[])std.zlib.uncompress(data, uncompressed_size, 31);
	    }
	    return uncompressed_data_;
	}
	return null; // TODO: implement bzip2 decompression
    }
}

struct FileHeader {
    private bio.sam.header.SamHeader _header;
    private string _headertext;

    this(CramBlock block) {
	enforce(block.content_type == CramBlock.ContentType.fileHeader);
	enforce(block.compression_method == CramBlock.CompressionMethod.raw);
	auto data = block.uncompressed_data;
	int len;
	read(data, &len);
	_headertext = cast(string)(data[0 .. len]);
    }

    bio.sam.header.SamHeader header() {
	if (_header is null) {
	    _header = new bio.sam.header.SamHeader(_headertext);
	    _headertext = null;
	}
	return _header;
    }
}

enum EncodingType : ubyte {
    none = 0,
    external = 1,
    golomb = 2,
    huffmanInt = 3,
    byteArrayLen = 4,
    byteArrayStop = 5,
    beta = 6,
    subExp = 7,
    golombRice = 8,
    gamma = 9
}

struct Encoding {
    EncodingType type;

    static struct ExternalParams { itf8 id; }
    static struct BetaParams {
	itf8 offset; itf8 length;
	itf8 read(ref BitStream bitstream) {
	    int result;
	    foreach (k; 0 .. length) {
		result <<= 1;
		if (bitstream.front) result += 1;
		bitstream.popFront();
	    }
	    return itf8(result - offset);
	}
    }
    static struct GammaParams {
	itf8 offset;
	itf8 read(ref BitStream bitstream) {
	    size_t n;
	    int result = 1;
	    while (!bitstream.front) {
		bitstream.popFront();
		++n;
	    }
	    bitstream.popFront();
	    while (n--) {
		result <<= 1;
		if (bitstream.front) result += 1;
		bitstream.popFront();
	    }
	    return itf8(result - offset);
	}
    }
    static struct GolombParams { itf8 offset; itf8 M; }
    static struct SubExponentialParams {
	itf8 offset; itf8 k;
	itf8 read(ref BitStream bitstream) {
	    int result;
	    size_t i;
	    while (bitstream.front) {
		++i; // count 1s
		bitstream.popFront();
	    }
	    bitstream.popFront(); // skip separator
	    size_t n_bits = i == 0 ? k : i + k - 1;
	    foreach (j; 0 .. n_bits) {
		result <<= 1;
		result += bitstream.front ? 1 : 0;
		bitstream.popFront();
	    }

	    result += 1 << (i + k - 1);
	    return itf8(result - offset);
	}
    }

    static struct HuffmanParams {
	itf8[] alphabet;
	itf8[] bitlengths;

	void buildTree() {
	    root = new Node;
	    auto data = alphabet.zip(bitlengths).array();
	    auto sorted = data.sort!((x, y) => x[1].value < y[1].value, SwapStrategy.stable);
	    int prev_len = -1;
	    ulong code = 0;
	    foreach (elem; sorted) {
		if (prev_len != -1) {
		    code += 1;
		    code <<= (elem[1] - prev_len);
		}
		//       // writeln("\t", elem[0], " => ", std.string.format("%0*b", elem[1].value, code));
		root.insert(elem[0], elem[1].value, code);
		prev_len = elem[1];
	    }
	}

	void printCodes() {
	    auto data = alphabet.zip(bitlengths).array();
	    auto sorted = data.sort!((x, y) => x[1].value < y[1].value, SwapStrategy.stable);
	    int prev_len = -1;
	    ulong code = 0;
	    foreach (elem; sorted) {
		if (prev_len != -1) {
		    code += 1;
		    code <<= (elem[1] - prev_len);
		}
	       // writeln("\t", elem[0], " => ", std.string.format("%0*b", elem[1].value, code));
		prev_len = elem[1];
	    }
	}

	private {
	    struct Node {
		Node*[2] children;
		itf8 value;

		bool is_leaf() const {
		    return children[0] is null && children[1] is null;
		}

		void insert(itf8 symbol, int bitlength, ulong code) {
		    Node* node = &this;
		    foreach (k; 0 .. bitlength) {
			int bit = !!(code & (1UL << (bitlength - k - 1)));
			auto next_node = node.children[bit];
			if (next_node is null) {
			    next_node = node.children[bit] = new Node;
			}
			node = next_node;
		    }
		    node.value = symbol;
		}
	    }

	    Node* root = null;
	}

	itf8 read(ref BitStream bitstream) {
	    auto node = root;
	    while (!node.is_leaf) {
		auto bit = bitstream.front;
		bitstream.popFront();
		node = node.children[bit];
	    }
	    return node.value;
	}
    }
    static struct ByteArrayLenParams { Encoding* len_encoding, val_encoding; }
    static struct ByteArrayStopParams { ubyte delim; itf8 external_id; }

    string toString() const {
	switch(type) {
	case EncodingType.none: return "None";
	case EncodingType.external:
	    return "External(#" ~ external.id.to!string ~ ")";
	case EncodingType.golomb, EncodingType.golombRice:
	    return "Golomb(offset=" ~ golomb.offset.to!string ~
		", M=" ~ golomb.M.to!string ~ ")";
	case EncodingType.beta:
	    return "Beta(offset=" ~ beta.offset.to!string ~
		", length=" ~ beta.length.to!string ~ ")";
	case EncodingType.gamma:
	    return "Gamma(offset=" ~ beta.offset.to!string ~ ")";
	case EncodingType.huffmanInt:
	    return "Huffman(alphabet=" ~ huffman.alphabet.to!string ~
		", bitlengths=" ~ huffman.bitlengths.to!string ~ ")";
	case EncodingType.subExp:
	    return "SubExponential(offset=" ~ sub_exp.offset.to!string ~
		", k=" ~ sub_exp.k.to!string ~ ")";
	case EncodingType.byteArrayLen:
	    return "ByteArrayLen(len_encoding=" ~ byte_array_len.len_encoding.toString() ~
		", val_encoding=" ~ byte_array_len.val_encoding.toString() ~ ")";
	case EncodingType.byteArrayStop:
	    return "ByteArrayStop(delim=" ~ byte_array_stop.delim.to!string ~
		", external_id=#" ~ byte_array_stop.external_id.to!string ~ ")";
	default:
	    return "NYI";
	}
    }

    union {
	ExternalParams external;
	GolombParams golomb;
	HuffmanParams huffman;
	ByteArrayLenParams byte_array_len;
	ByteArrayStopParams byte_array_stop;
	BetaParams beta;
	SubExponentialParams sub_exp;
	GammaParams gamma;
    }

    private bool readFromStream(S)(auto ref S stream) {
	itf8 n_bytes = void;
	if (!.read(stream, &n_bytes)) return false;
	switch(type) {
	case EncodingType.none: return true;
	case EncodingType.external: return .read(stream, &external.id);
	case EncodingType.golomb, EncodingType.golombRice:
	    return .read(stream, &golomb.offset) && .read(stream, &golomb.M);
	case EncodingType.huffmanInt:
	    auto result = .read(stream, huffman.alphabet) && .read(stream, huffman.bitlengths);
	    huffman.buildTree();
	    return result;
	case EncodingType.beta:
	    return .read(stream, &beta.offset) && .read(stream, &beta.length);
	case EncodingType.gamma:
	    return .read(stream, &gamma.offset);
	case EncodingType.subExp:
	    return .read(stream, &sub_exp.offset) && .read(stream, &sub_exp.k);
	case EncodingType.byteArrayLen:
	    byte_array_len.len_encoding = new Encoding;
	    byte_array_len.val_encoding = new Encoding;
	    return .read(stream, byte_array_len.len_encoding) &&
		.read(stream, byte_array_len.val_encoding);
	case EncodingType.byteArrayStop:
	    return .read(stream, &byte_array_stop.delim) &&
		.read(stream, &byte_array_stop.external_id);
	default:
	   // writeln("unknown encoding type, seeking ", n_bytes, " bytes");
	    return .seekCur(stream, n_bytes);
	}
    }
}

/// Reads information about encoding from a stream
bool read(S)(auto ref S stream, Encoding* encoding) {
    itf8 type = void;
    if (!read(stream, &type)) return false;
    encoding.type = cast(EncodingType)type.value;
    return encoding.readFromStream(stream);
}

struct BitStream {
    private {
	const(ubyte)[] _data;
	size_t _read_bytes;
	CramBlock[] _external_blocks;
	size_t[] _read_bytes_from_external;
	ubyte _curr_mask;
    }

    this(const(ubyte)[] data, CramBlock[] external_blocks) {
	_data = data;
	_curr_mask = 1 << 7;
	_read_bytes = 0;
	if (external_blocks !is null) {
	    _read_bytes_from_external.length = external_blocks.length;
	    _external_blocks = external_blocks;
	}
    }

    alias bool bit;

    bit front() @property const {
	return (_data[0] & _curr_mask) != 0;
    }

    bool empty() @property const {
	return _data.length == 0;
    }

    void popFront() {
	_curr_mask >>= 1;
	if (_curr_mask == 0) {
	    _data = _data[1 .. $];
	    _curr_mask = 1 << 7;
	    ++_read_bytes;
	}
    }

    size_t bytes_consumed() @property const {
	return _read_bytes;
    }

    itf8 read(Encoding encoding) {
	switch (encoding.type) {
	    case EncodingType.none:
		const(ubyte)[] bytes = _data[];
		itf8 result;
		.read(bytes, &result);
		_read_bytes += bytes.ptr - _data.ptr;
		_curr_mask = 1 << 7;
		return result;
	    case EncodingType.huffmanInt:
		return encoding.huffman.read(this);
	    case EncodingType.beta:
		return encoding.beta.read(this);
	    case EncodingType.gamma:
		return encoding.gamma.read(this);
	    case EncodingType.subExp:
		return encoding.sub_exp.read(this);
	    case EncodingType.external:
		auto id = encoding.external.id;
		auto bytes_read = _read_bytes_from_external[id];
		auto data = _external_blocks[id].uncompressed_data[bytes_read .. $];
		auto bitstream = BitStream(data, null);
		auto integer = bitstream.read(Encoding(EncodingType.none));
		auto n_bytes = bitstream.bytes_consumed;
		_read_bytes_from_external[id] += n_bytes;
		return integer;
	    default:
		stderr.writeln("TODO: implement ", encoding.type);
		return itf8.init;
	}
    }

    /// Note that these functions don't allocate but return a slice from uncompressed data
    const(ubyte)[] readArray(Encoding encoding, size_t len) {
	switch (encoding.type) {
	    case EncodingType.external:
		auto id = encoding.external.id;
		auto bytes_read = _read_bytes_from_external[id];
		_read_bytes_from_external[id] += len;
		return _external_blocks[id].uncompressed_data[bytes_read .. $][0 .. len];
	    default:
		stderr.writeln("TODO: implement ", encoding.type);
		return [];
	}
    }

    /// ditto
    const(ubyte)[] readArray(Encoding encoding) {
	switch (encoding.type) {
	    case EncodingType.byteArrayLen:
		auto length = this.read(*(encoding.byte_array_len.len_encoding));
		return this.readArray(*(encoding.byte_array_len.val_encoding), length);
	    case EncodingType.byteArrayStop:
		auto delim = encoding.byte_array_stop.delim;
		auto id = encoding.byte_array_stop.external_id;
		auto bytes_read = _read_bytes_from_external[id];
		auto start = bytes_read;
		while (_external_blocks[id].uncompressed_data[bytes_read] != delim)
		    ++bytes_read;
		auto result = _external_blocks[id].uncompressed_data[start .. bytes_read];
		_read_bytes_from_external[id] = bytes_read + 1;
		return result;
	    default:
		stderr.writeln("TODO: implement ", encoding.type);
		return [];
	}
    }
}

struct TagDefinition {
    private char[2] _name;
    char type;
    Encoding encoding;

    string name() @property const {
	return cast(string)_name[];
    }

    int id() @property const {
	int result = _name[0];
	result <<= 8;
	result |= _name[1];
	result <<= 8;
	result |= type;
	return result;
    }
}

final class CompressionHeader {
    private {
	static struct PreservationMap {
	    ubyte rn;
	    ubyte ap;
	    ubyte rr;
	    ubyte[5] sm;
	    ubyte[] td;
	    size_t[] starts;

	    char[4][5] bases;

	    char subst(char ref_base, ubyte code) {
		assert(code < 4);
		auto ref_base_code = Base5(ref_base).internal_code;
		return bases[ref_base_code][code];
	    }
	}

	PreservationMap _pm;

	Encoding[string] _encodings;
	Encoding[int] _tag_encodings;

	// TODO: can we assume fixed order here?
	bool readPreservationMap(S)(auto ref S stream) {
	    itf8 size_in_bytes = void;
	    itf8 n_elements = void;
	    if (!read(stream, &size_in_bytes) ||
		!read(stream, &n_elements)) return false;
	    ubyte[2] key = void;
	    foreach (i; 0 .. n_elements) {
		if (!read(stream, key[])) return false;
		switch (cast(string)key[]) {
		case "RN": read(stream, &_pm.rn); break;
		case "AP": read(stream, &_pm.ap); break;
		case "RR": read(stream, &_pm.rr); break;
		case "SM": read(stream, _pm.sm[]); break;
		case "TD": {
		    itf8 len = void;
		    if (!read(stream,  &len)) return false;
		    _pm.td.length = len;
		    if (!read(stream, _pm.td)) return false;
		    _pm.starts ~= 0;
		    foreach (size_t j, char c; _pm.td) {
			if (c == '\0')
			    _pm.starts ~= j + 1;
		    }
		    break;
		}
		default: assert(0);
		}
	    }

	    foreach (ref_base_code; 0 .. 5) {
		ubyte bc;
		foreach (k; 1 .. 5) {
		    if (bc == ref_base_code)
			++bc;
		    auto code = (_pm.sm[ref_base_code] >> (8 - 2 * k)) & 3;
		    _pm.bases[ref_base_code][code] = Base5.fromInternalCode(bc).asCharacter;
		    ++bc;
		}
	    }
	    return true;
	}

	// ditto
	bool readDataSeriesEncodings(S)(auto ref S stream) {
	    itf8 size_in_bytes = void;
	    itf8 n_elements = void;
	    if (!read(stream, &size_in_bytes) ||
		!read(stream, &n_elements)) return false;
	    ubyte[2] key = void;
	    Encoding encoding = void;
	    foreach (k; 0 .. n_elements) {
		if (!read(stream, key[])) return false;
		if (!read(stream, &encoding)) return false;
		_encodings[(cast(string)key[]).idup] = encoding;
	    }
	    return true;
	}

	bool readTagEncodings(S)(auto ref S stream) {
	    itf8 size_in_bytes = void;
	    itf8 n_elements = void;
	    if (!read(stream, &size_in_bytes) ||
		!read(stream, &n_elements)) return false;
	    itf8 key = void;
	    Encoding encoding = void;
	    foreach (k; 0 .. n_elements) {
		if (!read(stream, &key)) return false;
		if (!read(stream, &encoding)) return false;
		_tag_encodings[key] = encoding;
	    }
	    foreach (id; 0 .. _pm.starts.length - 1) {
		auto start = _pm.starts[id];
		auto stop = _pm.starts[id + 1] - 1;
		auto str = cast(string)_pm.td[start .. stop];
		_tag_definitions ~= new TagDefinition[str.length / 3];
		foreach (t; 0 .. str.length / 3) {
		    TagDefinition tag;
		    tag._name[] = str[t * 3 .. t * 3 + 2];
		    tag.type = str[t * 3 + 2];
		    tag.encoding = _tag_encodings[tag.id];
		    _tag_definitions.back[t] = tag;
		}
	    }
	    return true;
	}
    }

    CramRecordEncodingPack encodings;
    ReadFeatureEncodingPack read_feature_encodings;

    bool read_names_included() @property const { return !!_pm.rn; }
    bool AP_is_delta_encoded() @property const { return !!_pm.ap; }
    bool reference_required() @property const { return !!_pm.rr; }
    Encoding encoding(string key) const { return _encodings[key]; }

    Encoding tagEncoding(TagDefinition tag_definition) const {
	return _tag_encodings[tag_definition.id];
    }

    char getSubstitution(char ref_base, ubyte code) {
	return _pm.subst(ref_base, code);
    }

    private TagDefinition[][] _tag_definitions;

    const(TagDefinition)[] tags(size_t id) const {
	return _tag_definitions[id];
    }

    this(CramBlock block) {
	const(ubyte)[] data = block.uncompressed_data;
	bool ok = readPreservationMap(data); assert(ok);
	ok = readDataSeriesEncodings(data); assert(ok);
	ok = readTagEncodings(data); assert(ok);

	encodings = CramRecordEncodingPack(this);
	read_feature_encodings = ReadFeatureEncodingPack(this);
    }
}

final class MappedSliceHeader {
    private {
	itf8 _ref_id;
	itf8 _start_pos;
	itf8 _span;
	itf8 _n_records;
	itf8 _record_counter;
	itf8 _n_blocks;
	long[] _content_ids;
	itf8 _embed_ref_content_id;
	ubyte[16] _md5;
    }

    this(CramBlock block) {
	enforce(block.content_type == CramBlock.ContentType.mappedSliceHeader);
	enforce(block.compression_method == CramBlock.CompressionMethod.raw);
	auto data = block.uncompressed_data;
	bool ok = read(data, &_ref_id) &&
	    read(data, &_start_pos) &&
	    read(data, &_span) &&
	    read(data, &_n_records) &&
	    read(data, &_record_counter) &&
	    read(data, &_n_blocks);
	assert(ok);
	itf8 n_content_ids;
	itf8 content_id;
	read(data, &n_content_ids);
	_content_ids.length = n_content_ids;
	foreach (i; 0 .. n_content_ids) {
	    read(data, &content_id);
	    _content_ids[i] = content_id.value;
	}
	ok = read(data, &_embed_ref_content_id); assert(ok);
	ok = read(data, _md5[]); assert(ok);
    }

    long ref_id() const @property { return _ref_id; }
    long start_pos() const @property { return _start_pos; }
    long span() const @property { return _span; }
    long n_records() const @property { return _n_records; }
    long record_counter() const @property{ return _record_counter; }
    long n_blocks() const @property { return _n_blocks; }
    const(long)[] content_ids() const @property { return _content_ids; }
    long embed_ref_content_id() const @property { return _embed_ref_content_id; }
    ubyte[16] md5() const @property { return _md5; }
}

struct CramBlockRange {
    private {
	size_t _length;
	CramBlock _front;
	const(ubyte)[] _stream;

	this(const(ubyte[]) data, size_t n) {
	    _length = n;
	    _stream = data[];
	    prepareFront();
	}

	void prepareFront() {
	    with (_front) {
		bool ok =
		    read(_stream, cast(ubyte*)(&compression_method)) &&
		    read(_stream, cast(ubyte*)(&content_type)) &&
		    read(_stream, &content_id) &&
		    read(_stream, &compressed_size) &&
		    read(_stream, &uncompressed_size);
		assert(ok);

		// just slice instead of copying
		compressed_data = _stream[0 .. compressed_size];
		_stream = _stream[compressed_size .. $];
	    }
	}
    }

    size_t length() const { return _length; }

    alias opDollar = length;

    bool empty() const { return _length == 0; }

    ref CramBlock front() { return _front; }

    void popFront() {
	if (--_length > 0)
	    prepareFront();
    }
}

struct CramContainer {
    itf8 ref_id;
    itf8 start_pos;
    itf8 alignment_span;
    itf8 n_records;
    itf8 record_counter;
    itf8 n_bases;
    itf8 n_blocks;
    itf8[] landmarks;
    ubyte[] data;

    CramContainer dup() const {
	CramContainer c;
	c = cast()this;
	c.landmarks = landmarks.dup;
	c.data = data.dup;
	return c;
    }

    CramBlockRange blocks() const {
	return typeof(return)(data, n_blocks);
    }

    private CramSlice[] _slices;
    CramSlice[] slices() {
	if (_slices.length > 0)
	    return _slices;

	CompressionHeader compression;

	foreach(block; blocks) {
	    if (block.content_type == CramBlock.ContentType.compressionHeader) {
		compression = new CompressionHeader(block);
	    } else if (block.content_type == CramBlock.ContentType.mappedSliceHeader) {
		_slices.length += 1;
		_slices.back.compression = compression;
		_slices.back.header = new MappedSliceHeader(block);
	    } else if (block.content_type == CramBlock.ContentType.externalData) {
		_slices.back.external_blocks ~= block;
	    } else if (block.content_type == CramBlock.ContentType.coreData) {
		_slices.back.core_data = block.uncompressed_data;
	    }
	}
	return _slices;
    }
}

bool read(S)(auto ref S stream, out CramContainer container) {
    with(container) {
	int length;
	if (!read(stream, &length)) return false;
	data.length = length;

	return read(stream, &ref_id) &&
	    read(stream, &start_pos) &&
	    read(stream, &alignment_span) &&
	    read(stream, &n_records) &&
	    read(stream, &record_counter) &&
	    read(stream, &n_bases) &&
	    read(stream, &n_blocks) &&
	    read(stream, landmarks) &&
	    read(stream, data);
    }
}

struct CramContainerRange {
    private {
	std.stream.Stream _stream;
	CramContainer _front;
	bool _empty;

	this(std.stream.Stream stream) {
	    _stream = stream;
	    prepareFront();
	}

	void prepareFront() {
	    if (!read(_stream, _front))
		_empty = true;
	}
    }

    /// In case of parallel processing, .dup should be used
    /// on the front element.
    ref CramContainer front() { return _front; }

    bool empty() const { return _empty; }

    void popFront() { prepareFront(); }
}

class CramReader : IBamSamReader {
    private {
	string _filename;
	std.stream.Stream _stream;
	FastaReader _fasta;
	CramFileDefinition _def;
	std.parallelism.TaskPool _task_pool;

	alias Stream = std.stream.Stream;

	// stream with file definition stripped from the beginning
	Stream getContainerStream() {
	    if (_filename != null) {
		auto file = new std.stream.File(_filename);
		file.seekCur(CramFileDefinition.sizeof);
		return file;
	    } else {
		return _stream;
	    }
	}
    }

    this(string filename, string fasta_filename, std.parallelism.TaskPool task_pool=taskPool) {
	_filename = filename;
	_fasta = new FastaReader(fasta_filename);
	_task_pool = task_pool; // not used yet

	auto _file = new bio.core.utils.stream.File(filename);
	if (!read(_file, _def)) {
	    throw new CramReadError("Failed to read CRAM file definition");
	} else {
	    // TODO: check file magic
	    CramContainer first;
	    read(_file, first);
	    enforce(first.blocks.length == 1 &&
		    first.blocks.front.content_type == CramBlock.ContentType.fileHeader,
		    "First container must contain only SAM header");
	    _header = FileHeader(first.blocks.front).header;

	    // shameless copy-paste from bio.sam.reader
	    _reference_sequences = new ReferenceSequenceInfo[_header.sequences.length];
	    foreach (sq; _header.sequences) {
		auto seq = ReferenceSequenceInfo(sq.name, sq.length);
		auto n = cast(int)_reference_sequences.length;
		_reference_sequence_dict[sq.name] = n;
		_reference_sequences[_header.getSequenceIndex(seq.name)] = seq;
	    }
	}
    }

    CramContainerRange containers() {
	return typeof(return)(getContainerStream());
    }

    auto blocks() {
	return containers().map!(c => c.blocks()).joiner();
    }

    auto reads() {
	return containers().map!(c => cramRecords(c, _header)).joiner();
    }

    /// reads all containers; speed is limited by I/O performance
    size_t numberOfRecords() {
	return reduce!`a+b`(0UL, containers().map!(c => c.n_records));
    }

    // bunch of boilerplate below
    private SamHeader _header;
    bio.sam.header.SamHeader header() @property {
	return _header;
    }

    private {
	ReferenceSequenceInfo[] _reference_sequences;
	int[string] _reference_sequence_dict;
    }

    ///
    const(bio.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const {
	return _reference_sequences;
    }

    ///
    bool hasReference(string reference) {
	return null != (reference in _reference_sequence_dict);
    }

    ///
    bio.bam.reference.ReferenceSequence opIndex(string ref_name) {
	enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " is not present in the header");
	auto ref_id = _reference_sequence_dict[ref_name];
	return ReferenceSequence(null, ref_id, _reference_sequences[ref_id]);
    }

    string filename() @property const {
	return _filename;
    }

    /// no-op
    void assumeSequentialProcessing() {}

    ///
    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property {
	return inputRangeObject(bam_reads);
    }

    auto bam_reads() @property {
	return containers().map!(c => bamRecords(c, _header, this, _fasta)).joiner();
    }
}

struct CramBitFlags {
    itf8 flags;
    this(itf8 flags) { this.flags = flags; }
    this(int flags) { this.flags = itf8(flags); }
    bool is_paired()                     @property const { return cast(bool)(flags & 1); }
    bool is_proper_pair()                @property const { return cast(bool)(flags & 2); }
    bool is_unmapped()                   @property const { return cast(bool)(flags & 4); }
    bool mate_is_unmapped()              @property const { return cast(bool)(flags & 8); }
    bool is_reverse_strand()             @property const { return cast(bool)(flags & 0x10); }
    bool mate_is_reverse_strand()        @property const { return cast(bool)(flags & 0x20); }
    bool is_first_of_pair()              @property const { return cast(bool)(flags & 0x40); }
    bool is_second_of_pair()             @property const { return cast(bool)(flags & 0x80); }
    bool is_secondary_alignment()        @property const { return cast(bool)(flags & 0x100); }
    bool failed_quality_control()        @property const { return cast(bool)(flags & 0x200); }
    bool is_duplicate()                  @property const { return cast(bool)(flags & 0x400); }
}

struct CompressionBitFlags {
    itf8 flags;
    bool qualities_stored_as_array() @property const { return cast(bool)(flags & 1); }
    bool detached()                  @property const { return cast(bool)(flags & 2); }
    bool has_mate_downstream()       @property const { return cast(bool)(flags & 4); }
}

struct MateFlags {
    ubyte flags;
    this(itf8 flags) { this.flags = cast(ubyte)flags.value; }
    bool is_reverse_strand() @property const { return cast(bool)(flags & 1); }
    bool is_mapped()         @property const { return cast(bool)(flags & 2); }
}

struct EncodingPack(encodings...) {
    private {
	static string fields(encodings...)() {
	    string result;
	    foreach (enc; encodings) {
		result ~= "Encoding "~enc~";";
	    }
	    return result;
	}

	static string constructor(encodings...)() {
	    string result = "this(CompressionHeader compression) {";
	    foreach (enc; encodings) {
		result ~= enc~`=compression.encoding("`~enc~`");`;
	    }
	    result ~= "}";
	    return result;
	}
    }
    mixin(fields!encodings());
    mixin(constructor!encodings());
}

alias EncodingPack!("BA", "QS", "BS", "IN", "DL", "RS",
		    "SC", "PD", "HC", "FC", "FP") ReadFeatureEncodingPack;

alias EncodingPack!("BF", "CF", "RI", "RL", "AP", "RG", "RN", "MF", "NS",
		    "NP", "TS", "NF", "TL", "FN", "MQ", "BA", "QS") CramRecordEncodingPack;

struct ReadFeature {
    static enum Type : char {
	readBase = 'B',
	substitution = 'X',
	insertion = 'I',
	deletion = 'D',
	insertBase = 'i',
	qualityScore = 'Q',
	referenceSkip = 'N',
	softClip = 'S',
	padding = 'P',
	hardClip = 'H'
    }

    private {
	static string fields(specs...)() {
	    char[] result;
	    foreach (spec; specs)
		result ~= spec.Type.stringof ~ " " ~ spec.Name ~ ";";
	    return cast(string)result;
	}

	static string reader(specs...)() {
	    string result = "void read(ref BitStream bitstream, "
			    "ReadFeatureEncodingPack* encodings) {";
	    foreach (spec; specs) {
		result ~= spec.Name~`=cast(`~spec.Type.stringof ~`)`;
		if (isIntegral!(spec.Type) || is(spec.Type == char))
		    result ~= `bitstream.read(encodings.`~spec.Code~`);`;
		else
		    result ~= `bitstream.readArray(encodings.`~spec.Code~`);`;
	    }
	    result ~= "}";
	    return result;
	}
    }

    static struct F(string code, string name, T) {
	enum Code = code;
	enum Name = name;
	alias T Type;
    }

    static struct ReadFeatureImpl(specs...) {
	mixin(fields!specs());
	mixin(reader!specs());
    }

    alias ReadFeatureImpl!(F!("BA", "base", char), F!("QS", "quality", byte)) ReadBase;
    alias ReadFeatureImpl!(F!("BS", "code", byte)) Substitution;
    alias ReadFeatureImpl!(F!("IN", "bases", string)) Insertion;
    alias ReadFeatureImpl!(F!("DL", "length", int)) Deletion;
    alias ReadFeatureImpl!(F!("BA", "base", char)) InsertBase;
    alias ReadFeatureImpl!(F!("QS", "quality", byte)) QualityScore;
    alias ReadFeatureImpl!(F!("RS", "length", int)) ReferenceSkip;
    alias ReadFeatureImpl!(F!("SC", "bases", string)) SoftClip;
    alias ReadFeatureImpl!(F!("PD", "length", int)) Padding;
    alias ReadFeatureImpl!(F!("HC", "length", int)) HardClip;

    int position;
    Type type;
    union {
	ReadBase read_base;
	Substitution subst;
	Insertion insertion;
	Deletion deletion;
	InsertBase insert_base;
	QualityScore qual_score;
	ReferenceSkip ref_skip;
	SoftClip soft_clip;
	Padding padding;
	HardClip hard_clip;
    }

    string toString() {
	string result = "ReadFeature(position = " ~ position.to!string() ~
			", type = " ~ type.to!string ~ " -- ";
	final switch (type) {
	    case Type.readBase:      result ~= read_base.to!string();   break;
	    case Type.substitution:  result ~= subst.to!string();       break;
	    case Type.insertion:     result ~= read_base.to!string();   break;
	    case Type.deletion:      result ~= deletion.to!string();    break;
	    case Type.softClip:      result ~= soft_clip.to!string();   break;
	    case Type.hardClip:      result ~= hard_clip.to!string();   break;
	    case Type.padding:       result ~= padding.to!string();     break;
	    case Type.insertBase:    result ~= insert_base.to!string(); break;
	    case Type.qualityScore:  result ~= qual_score.to!string();  break;
	    case Type.referenceSkip: result ~= ref_skip.to!string();    break;
	}
	result ~= ")";
	return result;
    }

    static ReadFeature readFromBitStream(ref BitStream bitstream,
					 ReadFeatureEncodingPack* encodings,
					 ref int prev_pos)
    {
	ReadFeature feature;
	feature.type = cast(Type)(bitstream.read(encodings.FC));
	if (prev_pos == -1) {
	    // first position is 1-based;
	    feature.position = bitstream.read(encodings.FP) - 1;
	} else {
	    feature.position = bitstream.read(encodings.FP) + prev_pos;
	}
	prev_pos = feature.position;
	final switch (feature.type) {
	    case Type.readBase:      feature.read_base.read(bitstream, encodings);   break;
	    case Type.substitution:  feature.subst.read(bitstream, encodings);       break;
	    case Type.insertion:     feature.read_base.read(bitstream, encodings);   break;
	    case Type.deletion:      feature.deletion.read(bitstream, encodings);    break;
	    case Type.softClip:      feature.soft_clip.read(bitstream, encodings);   break;
	    case Type.hardClip:      feature.hard_clip.read(bitstream, encodings);   break;
	    case Type.padding:       feature.padding.read(bitstream, encodings);     break;
	    case Type.insertBase:    feature.insert_base.read(bitstream, encodings); break;
	    case Type.qualityScore:  feature.qual_score.read(bitstream, encodings);  break;
	    case Type.referenceSkip: feature.ref_skip.read(bitstream, encodings);    break;
	}
	return feature;
    }
}

struct CramRecord {
    CramBitFlags bit_flags;
    CompressionBitFlags compression_flags;
    int ref_id;
    int read_length;
    int position;
    int read_group_id;

    /// set if CRAM file includes all read names or if the mate is detached
    string read_name;

    /// set if compression_flags.detached is true
    MateFlags mate_flags;
    /// ditto
    int mate_ref_id;
    /// ditto
    int mate_position;
    /// ditto
    int template_size;

    /// set if compression_flags.detached is false
    int records_to_next_fragment;

    Tuple!(const(TagDefinition), const(ubyte)[])[] tags;

    /// set if bit_flags.is_unmapped is false
    ReadFeature[] features;
    /// ditto
    int mapping_quality;

    /// set if bit_flags.is_unmapped is true
    char[] bases;

    /// set if compression_flags.qualities_stored_as_array is true
    const(ubyte)[] qualities;
}

struct CramSlice {
    CompressionHeader compression;
    MappedSliceHeader header;
    const(ubyte)[] core_data;
    CramBlock[] external_blocks;
}

abstract class CramIterator {
    protected {
	CompressionHeader compression;
	MappedSliceHeader slice_header;
	CramBlock[] external_blocks;
	SamHeader _header;
    }

    private {
	const(ubyte)[] core_data;

	size_t current_record_index;
	size_t n_records;
	CramRecordEncodingPack* encodings;
	ReadFeatureEncodingPack* read_feature_encodings;
	BitStream bit_stream;
	int prev_pos;
    }

    this(CramSlice slice, SamHeader header) {
	_header = header;
	compression = slice.compression;
	slice_header = slice.header;
	core_data = slice.core_data;;
	external_blocks = slice.external_blocks;

	onInitializeHeaderBlocks();

	n_records = slice_header.n_records;
	if (n_records > 0) {
	    encodings = &compression.encodings;
	    read_feature_encodings = &compression.read_feature_encodings;

	    bit_stream = BitStream(core_data, external_blocks);
	    prev_pos = cast(int)slice_header.start_pos;

	    fetchNextRecord();
	}
    }

    bool empty() @property const {
	return n_records == current_record_index;
    }

    void popFront() {
	++current_record_index;
	if (current_record_index == n_records)
	    return;
	fetchNextRecord();
    }

    void onInitializeHeaderBlocks() {}

    abstract void setBitFlags(CramBitFlags bit_flags);
    abstract void setCompressionFlags(CompressionBitFlags compression_flags);
    abstract void setReferenceId(int ref_id);
    abstract void setReadLength(int read_length);
    abstract void setPosition(int position);
    abstract void setReadName(string read_name);
    abstract void setReadGroupId(int rg_id);
    abstract void setMateFlags(MateFlags mate_flags);
    abstract void setMateReferenceId(int ref_id);
    abstract void setMatePosition(int position);
    abstract void setTemplateSize(int template_size);
    abstract void setRecordsToNextFragment(int n_records);
    abstract void appendTag(const(TagDefinition) tag, const(ubyte)[] data);
    abstract void setNumberOfFeatures(size_t n_features);
    abstract void setFeature(size_t k, ReadFeature feature);
    abstract void setMappingQuality(int mapping_quality);
    abstract void setBase(size_t k, char base);
    abstract void setQualities(const(ubyte)[] qualities);
    abstract void finish();

    private void fetchNextRecord() {
	auto bit_flags = CramBitFlags(bit_stream.read(encodings.BF));
	setBitFlags(bit_flags);
	auto compression_flags = CompressionBitFlags(bit_stream.read(encodings.CF));
	setCompressionFlags(compression_flags);
	setReferenceId(bit_stream.read(encodings.RI));
	int read_length = bit_stream.read(encodings.RL);
	setReadLength(read_length);

	int position;
	if (!compression.AP_is_delta_encoded)
	    position = bit_stream.read(encodings.AP);
	else
	    position = prev_pos + bit_stream.read(encodings.AP);
	prev_pos = position;
	setPosition(position);

	setReadGroupId(bit_stream.read(encodings.RG));

	if (compression.read_names_included)
	    setReadName(cast(string)bit_stream.readArray(encodings.RN));

	if (compression_flags.detached) {
	    setMateFlags(MateFlags(bit_stream.read(encodings.MF)));
	    // logic behind this:
	    // even if we don't include read names in general, for detached mates
	    // they still should be present, since this is the way we find the mate later
	    if (!compression.read_names_included)
		setReadName(cast(string)bit_stream.readArray(encodings.RN));
	    setMateReferenceId(bit_stream.read(encodings.NS));
	    setMatePosition(bit_stream.read(encodings.NP));
	    setTemplateSize(bit_stream.read(encodings.TS));
	} else if (compression_flags.has_mate_downstream) {
	    setRecordsToNextFragment(bit_stream.read(encodings.NF));
	}

	int taglist_id = bit_stream.read(encodings.TL);
	foreach (ref tag; compression.tags(taglist_id)) {
	    appendTag(tag, bit_stream.readArray(tag.encoding));
	}

	if (!bit_flags.is_unmapped) {
	    int n_read_features = bit_stream.read(encodings.FN);
	    int feature_prev_pos = -1;
	    setNumberOfFeatures(n_read_features);
	    foreach (k; 0 .. n_read_features)
		setFeature(k, ReadFeature.readFromBitStream(bit_stream,
							    read_feature_encodings,
							    feature_prev_pos));
	    setMappingQuality(bit_stream.read(encodings.MQ));
	} else {
	    foreach (k; 0 .. read_length)
		setBase(k, cast(char)(bit_stream.read(encodings.BA)));
	    setMappingQuality(0xFF);
	}

	if (compression_flags.qualities_stored_as_array) {
	    setQualities(bit_stream.readArray(encodings.QS, read_length));
	} else {
	    setQualities([]);
	}

	finish();
    }
}

final class SliceCramRecordRange : CramIterator {
    private {
	CramRecord _front;
    }

    this(CramSlice slice, SamHeader header) {
	super(slice, header);
    }

    override void finish() {}

    override void popFront() {
	_front.tags.length = 0;
	super.popFront();
    }

    override void setBitFlags(CramBitFlags bit_flags) {
	_front.bit_flags = bit_flags;
    }

    override void setCompressionFlags(CompressionBitFlags compression_flags) {
	_front.compression_flags = compression_flags;
    }

    override void setReferenceId(int ref_id) {
	_front.ref_id = ref_id;
    }

    override void setReadLength(int read_length) {
	_front.read_length = read_length;
    }

    override void setPosition(int position) {
	_front.position = position;
    }

    override void setReadName(string read_name) {
	_front.read_name = read_name;
    }

    override void setReadGroupId(int rg_id) {
	_front.read_group_id = rg_id;
    }

    override void setMateFlags(MateFlags mate_flags) {
	_front.mate_flags = mate_flags;
    }

    override void setMateReferenceId(int ref_id) {
	_front.mate_ref_id = ref_id;
    }

    override void setMatePosition(int position) {
	_front.mate_position = position;
    }

    override void setTemplateSize(int template_size) {
	_front.template_size = template_size;
    }

    override void setRecordsToNextFragment(int n_records) {
	_front.records_to_next_fragment = n_records;
    }

    override void appendTag(const(TagDefinition) tag, const(ubyte)[] data) {
	_front.tags ~= tuple(tag, data);
    }

    override void setNumberOfFeatures(size_t n_features) {
	_front.features.length = n_features;
    }

    override void setFeature(size_t k, ReadFeature feature) {
	_front.features[k] = feature;
    }

    override void setMappingQuality(int mapping_quality) {
	_front.mapping_quality = mapping_quality;
    }

    override void setBase(size_t k, char base) {
	if (_front.bases.length != _front.read_length)
	    _front.bases.length = _front.read_length;
	_front.bases[k] = base;
    }

    override void setQualities(const(ubyte)[] qualities) {
	_front.qualities = qualities;
    }

    CramRecord front() @property {
	return _front;
    }
}

SliceCramRecordRange cramRecords(CramSlice slice, SamHeader header) {
    return new typeof(return)(slice, header);
}

auto cramRecords(ref CramContainer container, SamHeader header) {
    return container.slices.map!(s => cramRecords(s, header)).joiner();
}

void setBase(ubyte[] bases, size_t k, char base) {
    auto offset = k / 2;
    if (!(k & 1)) {
	bases[offset] &= 0x0F;
	bases[offset] = cast(ubyte)((Base16(base).internal_code << 4) | bases[offset]);
    } else {
	bases[offset] &= 0xF0;
	bases[offset] += cast(ubyte)(Base16(base).internal_code);
    }
}

void setBases(T)(ubyte[] dst, size_t offset, T[] src) {
    if (src.length == 0)
	return;
    if (offset % 2 != 0) {
	dst.setBase(offset, src[0]);
	src = src[1 .. $];
	++offset;
    }
    size_t real_offset = offset / 2;
    foreach (k; 0 .. src.length / 2) {
	auto b1 = Base16(src[2*k]).internal_code;
	auto b2 = Base16(src[2*k+1]).internal_code;
	dst[real_offset + k] = cast(ubyte)((b1 << 4) | b2);
    }
    if ((src.length & 1) != 0)
	dst[real_offset + src.length / 2] = cast(ubyte)(Base16(src[$-1]).internal_code << 4);
}

class FastaReader {
    private FastaSequence[string] _dict;

    this(string fn) {
	auto index_fn = fn ~ ".fai";

	auto fasta_file = new MmFile(fn);
	auto index = std.stdio.File(index_fn);
	foreach (line; index.byLine()) {
	    FastaSequence seq;
	    auto fields = splitter(line);
	    auto name = fields.front().idup; fields.popFront();
	    seq.length = fields.front().to!ulong; fields.popFront();
	    seq.offset = fields.front().to!ulong; fields.popFront();
	    seq.bases_per_line = fields.front().to!ulong; fields.popFront();
	    seq.bytes_per_line = fields.front().to!ulong; fields.popFront();
	    seq._file = fasta_file;
	    _dict[name] = seq;
	}
    }

    FastaSequence opIndex(string name) {
	return _dict[name];
    }
}

struct FastaSequence {
    ulong length;
    ulong offset;
    ulong bases_per_line;
    ulong bytes_per_line;

    private MmFile _file;

    private size_t byteNumber(size_t line_number, size_t inner_offset) {
	return offset + line_number * bytes_per_line + inner_offset;
    }

    char opIndex(size_t pos) {
	auto line_number = pos / bases_per_line;
	auto inner_offset = pos % bases_per_line;
	return cast(char)_file[byteNumber(line_number, inner_offset)];
    }

    void copySlice(size_t from, size_t to, ubyte* ptr) {
	auto first_line_number = from / bases_per_line;
	auto first_inner_offset = from % bases_per_line;
	auto first_byte_number = byteNumber(first_line_number, first_inner_offset);
	auto last_line_number = to / bases_per_line;
	auto last_inner_offset = to % bases_per_line;
	if (first_line_number == last_line_number) {
	    auto len = last_inner_offset - first_inner_offset;
	    ptr[0 .. len] = cast(ubyte[])_file[first_byte_number .. first_byte_number + len];
	} else {
	    auto len1 = bases_per_line - first_inner_offset;
	    ptr[0 .. len1] = cast(ubyte[])_file[first_byte_number .. first_byte_number + len1];
	    ptr += len1;
	    foreach (line_num; first_line_number + 1 .. last_line_number) {
		auto b = byteNumber(line_num, 0);
		ptr[0 .. bases_per_line] = cast(ubyte[])_file[b .. b + bases_per_line];
		ptr += bases_per_line;
	    }
	    auto len2 = last_inner_offset;
	    auto b = byteNumber(last_line_number, 0);
	    ptr[0 .. len2] = cast(ubyte[])_file[b .. b + len2];
	}
    }
}

struct ReadFeatureBuffer {
    private {
	ReadFeature[] _read_features_buf;
	ReadFeature[] _read_features() @property { return _read_features_buf[0 .. _curr]; }
	size_t _curr;
	ubyte[] _ref_chunk_buf;
    }

    FastaSequence reference;
    CompressionHeader* compression;

    void reserve(size_t n) {
	_read_features_buf.length = n;
	_curr = 0;
    }

    void append(ref ReadFeature feature) {
	assert(_curr < _read_features_buf.length);
	_read_features_buf[_curr++] = feature;
    }

    void reset() {
	_curr = 0;
    }

    ubyte[] restoreBases(ubyte[] bases, size_t position, size_t read_length) {
	int pos = 0;

	long ref_pos = position;
	long end_ref_pos = ref_pos + read_length;
	foreach (feature; _read_features) {
	    if (feature.type == 'I')
		end_ref_pos -= feature.insertion.bases.length;
	    else if (feature.type == 'i')
		end_ref_pos -= 1;
	    else if (feature.type == 'D')
		end_ref_pos += feature.deletion.length;
	    else if (feature.type == 'N')
		end_ref_pos += feature.ref_skip.length;
	    else if (feature.type == 'P')
		end_ref_pos += feature.padding.length;
	}

	end_ref_pos = min(reference.length, end_ref_pos);
	_ref_chunk_buf.length = max(_ref_chunk_buf.length, end_ref_pos - ref_pos);
	reference.copySlice(ref_pos, end_ref_pos, _ref_chunk_buf.ptr);
	auto ref_chunk = _ref_chunk_buf[0 .. end_ref_pos - ref_pos];

	if (_read_features.empty) {
	    if (reference.length < position + read_length) {
		auto overlap = reference.length - position;
		bases.setBases(0, ref_chunk[0 .. overlap]);
		foreach (k; overlap .. read_length)
		    bases.setBase(k, 'N');
	    } else {
		bases.setBases(0, ref_chunk[0 .. read_length]);
	    }
	}

	ref_pos = 0;

	foreach (feature; _read_features) {
	    bases.setBases(pos, ref_chunk[ref_pos .. ref_pos + feature.position - pos]);
	    pos = feature.position;

	    final switch (feature.type) {
	    case 'X':
		auto b = compression.getSubstitution(ref_chunk[ref_pos], feature.subst.code);
		bases.setBase(pos, b);
		++pos; ++ref_pos;
		break;
	    case 'I':
		bases.setBases(pos, feature.insertion.bases);
		pos += feature.insertion.bases.length;
		break;
	    case 'S':
		bases.setBases(pos, feature.soft_clip.bases);
		pos += feature.soft_clip.bases.length;
		break;
	    case 'H':
		break;
	    case 'N':
		ref_pos += feature.ref_skip.length;
		break;
	    case 'P':
		ref_pos += feature.padding.length;
		break;
	    case 'D':
		ref_pos += feature.deletion.length;
		break;
	    case 'i':
		bases.setBase(pos++, feature.insert_base.base);
		break;
	    case 'B':
		bases.setBase(pos++, feature.read_base.base);
		break;
	    case 'Q':
		break;
	    }
	}

	bases.setBases(pos, ref_chunk[ref_pos .. $]);

	return bases[0 .. read_length / 2 + read_length % 2];
    }

    void restoreQualityScores(ref ubyte[] scores) {
	foreach (feature; _read_features) {
	    switch (feature.type) {
	    case 'B':
		scores[feature.position] = feature.read_base.quality; break;
	    case 'Q':
		scores[feature.position] = feature.qual_score.quality; break;
	    default:
		break;
	    }
	}
    }

    CigarOperation[] restoreCigar(CigarOperation[] cigar, size_t read_length) {
	int last_op_len = 0;
	int last_op_pos = 0;
	int n_cigar_op = 0;
	char last_type = 'M';
	char co;

	foreach (feature; _read_features) {
	    int gap = feature.position - (last_op_pos + last_op_len);
	    if (gap > 0) {
		if (last_type != 'M') {
		    cigar[n_cigar_op++] = CigarOperation(last_op_len, last_type);
		    last_op_pos += last_op_len;
		    last_op_len = gap;
		} else {
		    last_op_len += gap;
		}

		last_type = 'M';
	    }

	    int len;

	    final switch (feature.type) {
	    case 'I': co = 'I'; len = cast(int)feature.insertion.bases.length; break;
	    case 'S': co = 'S'; len = cast(int)feature.soft_clip.bases.length; break;
	    case 'H': co = 'H'; len = feature.hard_clip.length; break;
	    case 'i': co = 'I'; len = 1; break;
	    case 'D': co = 'D'; len = feature.deletion.length; break;
	    case 'N': co = 'N'; len = feature.ref_skip.length; break;
	    case 'P': co = 'P'; len = feature.padding.length; break;
	    case 'X', 'B': co = 'M'; len = 1; break;
	    case 'Q': break;
	    }

	    if (last_type != co) {
		if (last_op_len > 0) {
		    cigar[n_cigar_op++] = CigarOperation(last_op_len, last_type);
		}
		last_type = co;
		last_op_len = len;
		last_op_pos = feature.position;
	    } else {
		last_op_len += len;
	    }

	    if (!CigarOperation(len, co).is_query_consuming)
		last_op_pos -= len;
	}

	if (last_type != 'M') {
	    cigar[n_cigar_op++] = CigarOperation(last_op_len, last_type);
	    if (read_length > last_op_pos + last_op_len) {
		auto len = read_length - (last_op_len + last_op_pos);
		cigar[n_cigar_op++] = CigarOperation(cast(uint)len, 'M');
	    }
	} else if (read_length > last_op_pos) {
	    auto len = read_length - last_op_pos;
	    cigar[n_cigar_op++] = CigarOperation(cast(uint)len, 'M');
	}

	if (n_cigar_op == 0) {
	    cigar[0] = CigarOperation(cast(uint)read_length, 'M');
	    return cigar[0 .. 1];
	}

	return cigar[0 .. n_cigar_op];
    }
}

final class SliceBamRecordRange : CramIterator {
    // based on net.sf.cram.encoding.reader.ReaderToBam   and
    //          net.sf.cram.encoding.reader.BAMRecordView from CramTools
    private {
	BamRead _front;
	CramReader _reader;
	FastaReader _fasta;

	uint recordCounter = 1;
	ubyte[] buf;
	int[] index;
	int[] distances;
	int[] names;
	int[] next;
	int[] prev;

	int flags;
	CompressionBitFlags compression_flags;
	int read_group_id;
	int mate_flags;

	int tag_data_len;
	ubyte[] tag_data;

	ubyte[] _ref;
	enum maxReadBufferLength = 1024 * 1024;
	ubyte[] bases;
	ubyte[] scores;
	ReadFeatureBuffer feature_buf;
	CigarOperation[] cigar_buf;
	size_t qual_len;

	const(ubyte)[][] read_groups;

	enum blockSizeOffset = 0;
	enum refIdOffset = 4;
	enum posOffset = 8;
	enum readNameLenOffset = 12;
	enum mappingQualityOffset = 13;
	enum indexBinOffset = 14;
	enum cigarLenOffset = 16;
	enum flagsOffset = 18;
	enum readLenOffset = 20;
	enum mateRefIdOffset = 24;
	enum matePosOffset = 28;
	enum insSizeOffset = 32;
	enum readNameOffset = 36;
	long cigarOffset = -1;
	long basesOffset = -1;
	long scoresOffset = -1;
	long endOffset = -1;

	long start;

	void reset() {
	    start = 0;
	    cigarOffset = -1;
	    basesOffset = -1;
	    scoresOffset = -1;
	    endOffset = -1;
	    tag_data_len = 0;
	    feature_buf.reset();
	}

	T read(T)(int at) const {
	    return *(cast(T*)(buf.ptr + start + at));
	}

	void write(T, string descr, U)(U value, int at) {
	    *cast(T*)(buf.ptr + start + at) = cast(T)value;
	    // writeln("[write ", descr, "] value = ", value, ", ", buf[start .. $][0 .. 32]);
	}
    }

    BamRead front() @property {
	return _front;
    }

    this(CramSlice slice, SamHeader header, CramReader reader,
	 FastaReader fasta) {
	this._reader = reader;
	this._fasta = fasta;
	buf = uninitializedArray!(ubyte[])(1024 * 1024 * 100);
	index.length = 4 * 100000;
	distances.length = 4 * 100000;
	names.length = 4 * 100000;
	next.length = distances.length;
	prev.length = distances.length;

	tag_data = uninitializedArray!(ubyte[])(1024 * 1024);

	bases = uninitializedArray!(ubyte[])(maxReadBufferLength);
	scores = uninitializedArray!(ubyte[])(maxReadBufferLength);
	cigar_buf = uninitializedArray!(CigarOperation[])(maxReadBufferLength);

	read_groups.length = header.read_groups.length;
	size_t i;
	foreach (rg; header.read_groups.values) {
	    read_groups[i] = cast(const(ubyte)[])rg.identifier;
	    while (read_groups[i].back == 0)
		read_groups[i] = read_groups[i][0 .. $ - 1];
	    ++i;
	}

	super(slice, header);
    }

    override void onInitializeHeaderBlocks() {
	if (slice_header.ref_id >= 0) {
	    auto seq_name = _reader.reference_sequences[slice_header.ref_id].name;
	    feature_buf.reference = _fasta[seq_name];
	}
	feature_buf.compression = &compression;
    }

    override void setBitFlags(CramBitFlags bit_flags) {
	this.flags = bit_flags.flags;
    }

    private void writeFlags() {
	write!(ushort, "flags")(flags, flagsOffset);
    }

    override void setCompressionFlags(CompressionBitFlags compression_flags) {
	this.compression_flags = compression_flags;
    }

    override void setReferenceId(int ref_id) { write!(int, "ref. id")(ref_id, refIdOffset); }

    override void setReadLength(int read_length) {
	write!(int, "read length")(read_length, readLenOffset);
    }

    private int getReadLength() const { return read!int(readLenOffset); }

    override void setPosition(int position) { write!(int, "position")(position-1, posOffset); }

    override void setReadName(string read_name) {
	write!(ubyte, "read name length")(read_name.length + 1, readNameLenOffset);
	auto from = start + readNameOffset;
	auto len = read_name.length;
	buf[from .. from + len] = cast(const(ubyte)[])read_name[];
	buf[from + len] = 0;
	cigarOffset = readNameOffset + len + 1;
    }

    bool read_name_is_set() @property const { return cigarOffset > -1; }

    void setIndexBin(int bin) { write!(ushort, "bin")(bin, indexBinOffset); }

    override void setReadGroupId(int rg_id) {
	read_group_id = rg_id;
    }

    // called only if mate is in the same slice
    override void setMateFlags(MateFlags mate_flags) {
	if (mate_flags.is_reverse_strand)
	    flags |= 0x20;
	if (!mate_flags.is_mapped)
	    flags |= 0x8;

	distances[recordCounter] = 0;
	next[recordCounter] = -1;
	prev[recordCounter] = -1;
    }

    override void setMateReferenceId(int ref_id) {
	write!(int, "mate ref. id")(ref_id, mateRefIdOffset);
    }

    override void setMatePosition(int position) {
	write!(int, "mate position")(position, matePosOffset);
    }

    override void setTemplateSize(int template_size) {
	write!(int, "template size")(template_size, insSizeOffset);
    }

    // called only if mate is not in the same slice
    override void setRecordsToNextFragment(int n_records) {
	mate_flags = 0;
	distances[recordCounter] = n_records;
	next[recordCounter] = recordCounter + n_records + 1; // <---.
	prev[next[recordCounter]] = recordCounter;         //        \
	names[recordCounter + n_records] = recordCounter; // TODO: ask Vadim about this
							  // inconsistency; is that a bug?
    }

    override void appendTag(const(TagDefinition) tag, const(ubyte)[] data) {
	tag_data[tag_data_len] = cast(ubyte)tag.name[0];
	tag_data[tag_data_len + 1] = cast(ubyte)tag.name[1];
	tag_data[tag_data_len + 2] = cast(ubyte)tag.type;
	tag_data_len += 3;

	bool is_string = tag.type == 'Z' || tag.type == 'H';
	if (is_string) {
	    while (!data.empty && data.back == 0)
		data = data[0 .. $ - 1];
	}

	tag_data[tag_data_len .. $][0 .. data.length] = data[];
	tag_data_len += data.length;

	if (is_string) {
	    tag_data[tag_data_len++] = 0;
	}
    }

    override void setNumberOfFeatures(size_t n_features) {
	feature_buf.reserve(n_features);
    }

    override void setFeature(size_t k, ReadFeature feature) {
	feature_buf.append(feature);
    }

    override void setMappingQuality(int mapping_quality) {
	write!(ubyte, "mapping quality")(mapping_quality, mappingQualityOffset);
    }

    override void setBase(size_t k, char base) {
	bases.setBase(k, base);
    }

    override void setQualities(const(ubyte)[] qualities) {
	if (qualities.length == 0) {
	    qual_len = getReadLength();
	    scores[0 .. qual_len] = 0xFF;
	} else {
	    qual_len = qualities.length;
	    scores[0 .. qual_len] = qualities[];
	}
    }

    override void finish() {
	if (!read_name_is_set) {
	    size_t int_name;
	    if (names[recordCounter] == 0) {
		int_name = recordCounter;
	    } else {
		int_name = names[recordCounter];
	    }
	    char[16] tmp = void;
	    char* ptr = tmp.ptr;
	    bio.core.utils.format.write(ptr, int_name);
	    setReadName(cast(string)(tmp[0 .. ptr - tmp.ptr]));
	}

	if (read_group_id >= 0) {
	    auto rg = cast(ubyte[])(read_groups[read_group_id]);
	    ubyte* ptr = tag_data.ptr + tag_data_len;
	    *ptr++ = 'R'; *ptr++ = 'G'; *ptr++ = 'Z';
	    ptr[0 .. rg.length] = rg[];
	    ptr[rg.length] = 0;
	    tag_data_len += 4 + rg.length;
	}

	auto read_length = getReadLength();

	// restore bases and CIGAR
	CigarOperation[] cigar;
	ubyte[] read_bases;
	if (!CramBitFlags(flags).is_unmapped) {
	    auto position = read!int(posOffset);
	    read_bases = feature_buf.restoreBases(bases, position, read_length);
	    cigar = feature_buf.restoreCigar(cigar_buf, read_length);
	    feature_buf.restoreQualityScores(scores);
	} else {
	    auto bases_len = read_length / 2 + read_length % 2;
	    read_bases = bases[0 .. bases_len];
	}

	// copy CIGAR
	basesOffset = cigarOffset + cigar.length * CigarOperation.sizeof;
	foreach (p; 0 .. cigar.length) {
	    auto offset = start + cigarOffset + p * CigarOperation.sizeof;
	    auto ptr = cast(CigarOperation*)(buf.ptr + offset);
	    *ptr = cigar[p];
	}

	// set CIGAR length
	write!(ushort, "cigar length")(cigar.length, cigarLenOffset);

	// copy bases
	buf[start + basesOffset .. $][0 .. read_bases.length] = read_bases[];

	// copy qualities
	scoresOffset = basesOffset + read_length / 2 + read_length % 2;
	buf[start + scoresOffset .. $][0 .. qual_len] = scores[0 .. qual_len];

	// copy tag data
	auto tagsOffset = scoresOffset + qual_len;
	endOffset = tagsOffset;
	// TODO: exception handling in case of overflow; use outbuffer?
	// writeln("[write] tags offset = ", tagsOffset - 4);
	// writeln("[write] qual offset = ", scoresOffset - 4);
	// writeln("[write] cigar offset = ", cigarOffset - 4);
	// writeln("[write] seq offset = ", basesOffset - 4);
	// writeln(tag_data[0 .. tag_data_len]);
	buf[start + tagsOffset .. $][0 .. tag_data_len] = tag_data[0 .. tag_data_len];
	endOffset += tag_data_len;

	// FIXME: this is temporary
	setMateReferenceId(-1);

	writeFlags();

	_front.raw_data = buf[start + refIdOffset .. start + endOffset];
	_front.associateWithReader(_reader);

	auto new_start = start + endOffset;
	reset();
	start = new_start;

	recordCounter += 1;
    }
}

auto bamRecords(CramSlice slice, SamHeader header,
		     CramReader reader, FastaReader fasta)
{
    return new SliceBamRecordRange(slice, header, reader, fasta);
}

auto bamRecords(ref CramContainer container, SamHeader header,
		CramReader reader, FastaReader fasta) {
    return container.slices.map!(s => bamRecords(s, header, reader, fasta)).joiner();
}

void main(string[] args) {
    import std.datetime;
    StopWatch sw;
    sw.start();
    auto cram = new CramReader(args[1], args[2]);
    size_t records_read;
    // scope(exit) stderr.writeln(records_read);
    int asdf = 1;
    auto writer = stdout.lockingTextWriter;
    import std.format;
    FormatSpec!char fmt;
    fmt.spec = 's';
    foreach (record; cram.bam_reads) {
	++records_read;
	record.toString((const(char)[] x) { writer.put(x); writer.put("\n"); }, fmt);
	if (records_read >= asdf * 1_000_000) {
	    stderr.writeln("Records read: ", records_read, ", time: ", sw.peek().msecs, " ms");
	    ++asdf;
	}
    }
}
