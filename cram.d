import std.stdio, std.stream, std.bitmanip, std.traits, std.exception, std.conv,
    std.algorithm, std.array, std.string, std.zlib;
import bio.sam.header;

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
    static struct BetaParams { itf8 offset; itf8 length; }
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
    static struct SubExponentialParams { itf8 offset; itf8 k; }

    static struct HuffmanParams {
        itf8[] alphabet; 
        union {
            itf8[] bitlengths;
            itf8[] codes;
        }

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
                //        writeln("\t", elem[0], " => ", std.string.format("%0*b", elem[1].value, code));
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
                writeln("\t", elem[0], " => ", std.string.format("%0*b", elem[1].value, code));
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
                enforce(!bitstream.empty, "error reading huffman code");
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
            writeln("unknown encoding type, seeking ", n_bytes, " bytes");
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
        ubyte _curr_byte_read_bits;
        CramBlock[] _external_blocks;
        size_t[] _read_bytes_from_external;
    }

    this(const(ubyte)[] data, CramBlock[] external_blocks) {
        _data = data;
        _curr_byte_read_bits = 0;
        _read_bytes_from_external.length = external_blocks.length;
        _external_blocks = external_blocks;
    }

    alias bool bit;

    bit front() @property const {
        return (_data[0] & (1 << (7 - _curr_byte_read_bits))) != 0;
    }

    bool empty() @property const {
        return _data.length == 0;
    }

    void popFront() {
        ++_curr_byte_read_bits;
        if (_curr_byte_read_bits == 8) {
            _data = _data[1 .. $];
            _curr_byte_read_bits = 0;
        }
    }

    itf8 read(Encoding encoding) {
        switch (encoding.type) {
            case EncodingType.huffmanInt:
                return encoding.huffman.read(this);
            case EncodingType.gamma:
                return encoding.gamma.read(this);
            default:
                stderr.writeln("TODO: implement ", encoding.type);
                return itf8.init;
        }
    }

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

    const(ubyte)[] readArray(Encoding encoding) {
        switch (encoding.type) {
            case EncodingType.byteArrayLen:
                auto length = this.read(*(encoding.byte_array_len.len_encoding));
                return this.readArray(*(encoding.byte_array_len.val_encoding), length);
            default:
                stderr.writeln("TODO: implement ", encoding.type);
                return [];
        }
    }
}

struct TagDefinition {
    private char[2] _name;
    char type;

    this(string name, char type) {
        _name[0] = name[0];
        _name[1] = name[1];
        this.type = type;
    }

    string name() @property const {
        return cast(string)_name[];
    }
}

struct CompressionHeader {
    static struct TagDefinitionRange {
        string s;
        this(string s) {
            this.s = s;
        }

        bool empty() @property const { return s.length == 0; }
        TagDefinition front() @property const {
            return TagDefinition(s[0 .. 2], s[2]);
        }
        void popFront() { s = s[3 .. $]; }
    }

    private {
        static struct PreservationMap {
            ubyte rn;
            ubyte ap;
            ubyte rr;
            ubyte[5] sm;
            ubyte[] td;
            size_t[] starts;
        }

        static struct SubstitutionMatrix {
            ubyte[5] sm;
            this(ubyte[5] sm) { this.sm = sm; }
        }

        PreservationMap _pm;

        Encoding[string] _encodings;

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
    }

    bool read_names_included() @property const { return !!_pm.rn; }
    bool AP_is_delta_encoded() @property const { return !!_pm.ap; }
    bool reference_required() @property const { return !!_pm.rr; }
    Encoding encoding(string key) const { return _encodings[key]; }

    auto tags(size_t id) const {
        auto start = _pm.starts[id];
        auto stop = _pm.starts[id + 1] - 1;
        return TagDefinitionRange(cast(string)_pm.td[start .. stop]);
    }
    
    this(CramBlock block) {
        const(ubyte)[] data = block.uncompressed_data;
        bool ok = readPreservationMap(data); assert(ok);
        ok = readDataSeriesEncodings(data); assert(ok);
    }
}

struct MappedSliceHeader {
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

class CramReader {
    private {
        string _filename;
        std.stream.Stream _stream;
        CramFileDefinition _def;

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
    
    this(string filename) {
        _filename = filename;

        auto _file = new std.stream.File(filename);
        if (!read(_file, _def)) {
            throw new CramReadError("Failed to read CRAM file definition");
        } else {
            // TODO: check file magic
        }
    }

    CramContainerRange containers() {
        return typeof(return)(getContainerStream());
    }

    auto blocks() {
        return containers().map!(c => c.blocks()).joiner();
    }
}

struct CramBitFlags {
    itf8 flags;
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

import std.range, std.algorithm;

void main(string[] args) {
    auto cram = new CramReader(args[1]);
    foreach (container; cram.containers.take(2)) {
        writeln("=== container ===", " (", container.data.length, " bytes, ", container.n_records, " records)");
        CompressionHeader compression;
        MappedSliceHeader slice_header;
        CramBlock[] external_blocks;
        const(ubyte)[] core_data;
        size_t external_id = 0;
        foreach(block; container.blocks) {
            if (block.content_type == CramBlock.ContentType.compressionHeader) {
                compression = CompressionHeader(block);
                writeln("BF: ", compression.encoding("BF"));
                compression.encoding("BF").huffman.printCodes();
                writeln("CF: ", compression.encoding("CF"));
                compression.encoding("CF").huffman.printCodes();
                writeln("RI: ", compression.encoding("RI"));
                compression.encoding("RI").huffman.printCodes();
                writeln("RL: ", compression.encoding("RL"));
                compression.encoding("RL").huffman.printCodes();
                writeln("AP: ", compression.encoding("AP"));
                writeln("RG: ", compression.encoding("RG"));
                writeln("QS: ", compression.encoding("QS"));
                writeln("RN: ", compression.encoding("RN"));
            } else if (block.content_type == CramBlock.ContentType.mappedSliceHeader) {
                slice_header = MappedSliceHeader(block);
            } else if (block.content_type == CramBlock.ContentType.externalData) {
                writeln("External block #", external_id, ": ", block.compression_method);
                ++external_id;
                writeln("Data size: ", block.uncompressed_size);
                writeln("Data: ", cast(string)block.uncompressed_data[0 .. 80], "...");
                external_blocks ~= block;
            } else if (block.content_type == CramBlock.ContentType.coreData) {
                core_data = block.uncompressed_data;
            }
        }

        if (container.n_records > 0) {
            auto bit_stream = BitStream(core_data, external_blocks);
            for (size_t i = 0; i < min(1, container.n_records); ++i) {
            auto bit_flags = CramBitFlags(bit_stream.read(compression.encoding("BF")));
            auto compression_flags = CompressionBitFlags(bit_stream.read(compression.encoding("CF")));
            int ref_id = bit_stream.read(compression.encoding("RI"));
            int read_length = bit_stream.read(compression.encoding("RL"));

            // compression.AP_is_delta_encoded?
            int position = bit_stream.read(compression.encoding("AP"));

            int read_group = bit_stream.read(compression.encoding("RG"));

            string read_name;
            if (compression.read_names_included)
                read_name = cast(string)bit_stream.readArray(compression.encoding("RN"));
            writeln("Read: ", read_name);
            writeln("  Position: ", position + container.start_pos);
            writeln("  Read group id: ", read_group);

            if (compression_flags.detached) {
                auto mate_flags = MateFlags(bit_stream.read(compression.encoding("MF")));
                string mate_read_name;
                if (compression.read_names_included)
                    mate_read_name = cast(string)bit_stream.readArray(compression.encoding("RN"));
                int mate_ref_id = bit_stream.read(compression.encoding("NS"));
                int mate_position = bit_stream.read(compression.encoding("NP"));
                int template_size = bit_stream.read(compression.encoding("TS"));
                writeln("Mate: ", mate_read_name);
                writeln("  Ref. id: ", mate_ref_id);
                writeln("  Position: ", mate_position);
                writeln("  Template size: ", template_size);
            } else if (compression_flags.has_mate_downstream) {
                int records_to_next_fragment = bit_stream.read(compression.encoding("NF"));
                writeln("Records to next fragment: ", records_to_next_fragment);
            }

            int taglist_id = bit_stream.read(compression.encoding("TL"));
            writeln("Tags: ", compression.tags(taglist_id));
        }
        }
    }
}

