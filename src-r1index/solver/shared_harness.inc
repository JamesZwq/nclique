// Frozen input, width and hash helpers.
static cpp_int binomial(unsigned n,unsigned k) {
    if(k>n)return 0;
    k=std::min(k,n-k);cpp_int value=1;
    for(unsigned i=1;i<=k;++i){value*=n-i+1;value/=i;}
    return value;
}

struct Input {
    Graph graph;
    Vertices ordinary;
    Vertex d=0;
    double load_ms=0,order_ms=0;
};

static Input prepare(const std::string& path) {
    Input result;auto start=Clock::now();auto raw=Graph::read(path);result.load_ms=ms(start);
    start=Clock::now();Seeds seeds(raw);result.d=seeds.maximum;result.ordinary.resize(raw.n);
    std::vector<std::pair<Vertex,Vertex>> edges;edges.reserve(raw.m);
    for(Vertex v=0;v<raw.n;++v) {
        result.ordinary[seeds.rank[v]]=seeds.ordinary[v];
        for(Vertex u:raw.row(v))if(v<u)edges.emplace_back(seeds.rank[v],seeds.rank[u]);
    }
    result.graph=Graph::from_edges(raw.n,std::move(edges));result.order_ms=ms(start);return result;
}

static cpp_int count_bound(const Graph& graph,const Layout& index,Vertex d) {
    std::vector<cpp_int> peaks(d+2),bound(graph.n);
    for(Vertex q=0;q<=d+1;++q)peaks[q]=binomial(q,q/2);
    for(size_t p=0;p<index.paths.size();++p) {
        const size_t q=index.paths.row(p).size()-index.paths.holds[p];
        require(q<peaks.size(),"path exceeds ordinary-core bound");
        for(Vertex v:index.paths.row(p))bound[v]+=peaks[q];
    }
    cpp_int maximum=d;
    for(const auto& value:bound)if(value>maximum)maximum=value;
    return maximum;
}

static unsigned width(const cpp_int& bound) {
    for(unsigned bits:{64,128,256,512})if(bound<(cpp_int(1)<<bits)-1)return bits;
    throw std::overflow_error("full-range bound exceeds 512-bit kernel");
}

template<class T> static uint64_t hash_core(const std::vector<T>& values) {
    uint64_t hash=1469598103934665603ULL;
    const T mask=std::numeric_limits<uint64_t>::max();
    for(T value:values) {
        for(unsigned i=0;i<std::numeric_limits<T>::digits;i+=64) {
            hash^=static_cast<uint64_t>(value&mask);hash*=1099511628211ULL;
            if constexpr(std::numeric_limits<T>::digits>64)value>>=64;
        }
    }
    return hash;
}

