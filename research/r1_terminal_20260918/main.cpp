#include "terminal.hpp"
#include "../r1_floorbatch_20260916/oracle.hpp"
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>

using namespace orderreplay;
#include "../r1_orderreplay_20260917/shared_harness.inc"

static Graph complete(Vertex n) {
    std::vector<std::pair<Vertex,Vertex>> edges;
    for(Vertex u=0;u<n;++u)for(Vertex v=u+1;v<n;++v)edges.emplace_back(u,v);
    return Graph::from_edges(n,std::move(edges));
}

static Graph split_graph(Vertex common,Vertex choices,bool extra_edge=false) {
    const Vertex n=common+choices+(extra_edge?2:0);
    std::vector<std::pair<Vertex,Vertex>> edges;
    for(Vertex u=0;u<common;++u)for(Vertex v=u+1;v<n;++v)edges.emplace_back(u,v);
    if(extra_edge)edges.emplace_back(n-2,n-1);
    return Graph::from_edges(n,std::move(edges));
}

static void coverage_and_residuals(const Graph& graph,const terminal::Index& index,uint64_t& checks) {
    const Vertex n=graph.n;
    for(int s=2;s<=index.maximum;++s) {
        auto wanted=bottomup::clique_masks(graph,s);
        std::sort(wanted.begin(),wanted.end());
        std::vector<uint64_t> covered;
        for(const auto& row:index.rows)if(row.valid(s)) {
            uint64_t hm=0,qm=0,xm=0;
            for(terminal::Offset i=row.begin;i<row.end;++i) {
                auto& mask=i<row.hold_end?hm:(i<row.pivot_end?qm:xm);
                mask|=uint64_t{1}<<index.members[i];
            }
            require(std::popcount(hm|qm|xm)==static_cast<int>(row.end-row.begin),"repeated terminal member");
            for(uint64_t mask=0;mask<(uint64_t{1}<<n);++mask) {
                if(std::popcount(mask)!=s || (mask&hm)!=hm || (mask&~(hm|qm|xm)))continue;
                const int choices=std::popcount(mask&xm);
                if(row.group==absent ? choices!=0 : (index.zero_choice[row.group]?choices>1:choices!=1))continue;
                covered.push_back(mask);
            }
        }
        std::sort(covered.begin(),covered.end());
        require(covered==wanted,"terminal clique coverage or multiplicity mismatch");++checks;
        if(n>5)continue;
        terminal::Solver<uint64_t>::Choose choose(n,index.maximum);
        for(uint64_t live=0;live<(uint64_t{1}<<n);++live) {
            std::vector<uint64_t> actual(n),count(n);
            for(uint64_t mask:wanted)if((mask&live)==mask)
                for(Vertex v=0;v<n;++v)if((mask>>v)&1)++actual[v];
            terminal::Metrics metrics;
            for(terminal::RowId p=0;p<index.rows.size();++p) {
                const auto& row=index.rows[p];if(!row.valid(s))continue;
                bool holds=true;Vertex q=0,z=0;
                for(terminal::Offset i=row.begin;i<row.end;++i) {
                    const bool inside=(live>>index.members[i])&1;
                    if(i<row.hold_end)holds&=inside;
                    else if(i<row.pivot_end)q+=inside;else z+=inside;
                }
                if(!holds)continue;
                const auto value=terminal::Solver<uint64_t>::coefficients(index,p,s,q,z,choose,metrics);
                for(terminal::Offset i=row.begin;i<row.end;++i)if((live>>index.members[i])&1)
                    count[index.members[i]]+=i<row.hold_end?value.h:(i<row.pivot_end?value.q:value.x);
            }
            require(count==actual,"arbitrary residual count mismatch");++checks;
        }
    }
}

static void check_replay(const Graph& graph,const Layout& original,const terminal::Index& index,
                         const std::vector<uint64_t>& truth,uint64_t& checks) {
    Kernel<uint64_t>::Combinations choose(graph.n,index.maximum);
    Vertices order(graph.n);std::iota(order.begin(),order.end(),0);
    std::mt19937_64 random(2026091822+graph.n+graph.m);
    for(int trial=0;trial<2;++trial) {
        if(trial)std::shuffle(order.begin(),order.end(),random);
        Vertices rank(graph.n);for(Vertex i=0;i<graph.n;++i)rank[order[i]]=i;
        orderreplay::Metrics metrics;
        terminal::ReplayPlan<uint64_t> plan(graph.n,index,choose,true,metrics);
        for(int s=2;s<=index.maximum;++s) {
            const auto cliques=bottomup::clique_masks(graph,s);
            std::vector<uint64_t> forward(graph.n),upper(graph.n),actual(graph.n),certificate(graph.n);
            for(auto mask:cliques) {
                Vertex first=absent;
                for(Vertex v=0;v<graph.n;++v)if((mask>>v)&1)
                    if(first==absent || rank[v]<rank[first])first=v;
                ++forward[first];
            }
            uint64_t level=0;
            for(Vertex v:order)if(original.initial_top[v]>=static_cast<Vertex>(s)) {
                level=std::max(level,forward[v]);upper[v]=level;
            }
            plan.forward(s,order,actual);
            require(actual==upper && plan.counts()==forward,"factored compiled forward mismatch");
            for(auto mask:cliques)for(Vertex v=0;v<graph.n;++v)if((mask>>v)&1) {
                bool inside=true;
                for(Vertex u=0;u<graph.n;++u)if(((mask>>u)&1) && upper[u]<upper[v])inside=false;
                certificate[v]+=inside;
            }
            const bool accepted=plan.trial(s,order,actual);
            require(actual==upper && plan.counts()==certificate,"factored compiled threshold mismatch");
            bool exact=true;
            for(Vertex v=0;v<graph.n;++v)exact&=upper[v]==truth[static_cast<size_t>(s)*graph.n+v];
            require(accepted==exact,"factored replay decision mismatch");++checks;
        }
    }
}

template<class T> static void check_small(const Graph& graph,int maximum,bool coverage,uint64_t& cells,uint64_t& families) {
    Seeds seeds(graph);const auto expected64=floorbatch::explicit_core(graph,maximum,graph.n<=5);
    const std::vector<T> expected(expected64.begin(),expected64.end());
    Layout original(graph,maximum);original.prepare(graph.n);
    typename Kernel<T>::Combinations choose(graph.n,maximum);
    Kernel<T> baseline;
    const auto fixed=baseline.template fixed_sparse<true>(original,graph.n,choose,seeds.ordinary);
    require(fixed.core==expected,"original reference differs from explicit oracle");
    const auto reference=baseline.template solve<true>(graph,original,choose,seeds.ordinary);
    require(reference.common.data.core==expected,"stream differs from explicit oracle");
    for(int mode=0;mode<3;++mode) {
        terminal::Index index(maximum);terminal::build(graph,index,mode);index.prepare(graph.n);
        if(coverage)coverage_and_residuals(graph,index,families);
        if(mode==0)require(index.members.size()==original.paths.vertices.size() && index.rows.size()==original.paths.size(),"plain direct completion changed paths");
        else require(index.members.size()+index.work.saved_members==original.paths.vertices.size(),"membership reduction identity differs");
        const auto result=terminal::Solver<T>::template solve<true>(graph,index,choose,seeds.ordinary,&expected);
        require(result.common.data.core==expected,"factored final core differs");cells+=expected.size();
        require(result.common.data.work.events==reference.common.data.work.events,"factoring changed source event count");
        const auto replay=terminal::Solver<T>::template solve<true,true>(graph,index,choose,seeds.ordinary,&expected);
        require(replay.common.data.core==expected,"factored replay final core differs");cells+=expected.size();
        if constexpr(std::is_same_v<T,uint64_t>)check_replay(graph,original,index,expected64,families);
    }
}

template<class T> static void check_wide(const Graph& graph,bool analytic) {
    Seeds seeds(graph);const int maximum=seeds.maximum+1;
    Layout original(graph,maximum);original.prepare(graph.n);
    typename Kernel<T>::Combinations choose(graph.n,maximum);
    const auto expected=Kernel<T>{}.template fixed_sparse<true>(original,graph.n,choose,seeds.ordinary).core;
    for(int mode=0;mode<3;++mode) {
        terminal::Index index(maximum);terminal::build(graph,index,mode);index.prepare(graph.n);
        const auto result=terminal::Solver<T>::solve(graph,index,choose,seeds.ordinary);
        require(result.common.data.core==expected,"wide factored output mismatch");
        const auto replay=terminal::Solver<T>::template solve<false,true>(graph,index,choose,seeds.ordinary);
        require(replay.common.data.core==expected,"wide factored replay mismatch");
    }
    if(analytic)for(int s=2;s<=maximum;++s)for(Vertex v=0;v<graph.n;++v)
        require(cpp_int(expected[static_cast<size_t>(s)*graph.n+v])==binomial(graph.n-1,s-1),"analytic complete-graph answer differs");
}

static void selftest() {
    uint64_t graphs=0,cells=0,families=0;std::mt19937_64 random(2026091821);
    for(Vertex n=0;n<=6;++n) {
        std::vector<std::pair<Vertex,Vertex>> pairs;
        for(Vertex u=0;u<n;++u)for(Vertex v=u+1;v<n;++v)pairs.emplace_back(u,v);
        for(uint64_t mask=0;mask<(uint64_t{1}<<pairs.size());++mask) {
            std::vector<std::pair<Vertex,Vertex>> edges;
            for(size_t i=0;i<pairs.size();++i)if((mask>>i)&1)edges.push_back(pairs[i]);
            const auto graph=Graph::from_edges(n,std::move(edges));Seeds seeds(graph);
            const int maximum=std::max(2,static_cast<int>(seeds.maximum)+1);
            check_small<uint64_t>(graph,maximum,true,cells,families);++graphs;
            if(n<=4 || mask%257==0) {
                check_small<unsigned __int128>(graph,maximum,false,cells,families);
                check_small<boost::multiprecision::uint256_t>(graph,maximum,false,cells,families);
            }
        }
    }
    for(int trial=0;trial<180;++trial) {
        const Vertex n=7+random()%6;const unsigned density=random()%101;
        std::vector<std::pair<Vertex,Vertex>> edges;
        for(Vertex u=0;u<n;++u)for(Vertex v=u+1;v<n;++v)if(random()%100<density)edges.emplace_back(u,v);
        const auto graph=Graph::from_edges(n,std::move(edges));Seeds seeds(graph);
        const int maximum=trial%3==0?16:std::max(2,static_cast<int>(seeds.maximum)+1);
        check_small<uint64_t>(graph,maximum,true,cells,families);++graphs;
    }
    for(bool extra:{false,true})for(Vertex h:{1,2,4}) {
        const auto graph=split_graph(h,4,extra);
        check_small<uint64_t>(graph,graph.n,true,cells,families);++graphs;
    }
    check_wide<uint64_t>(complete(34),true);
    check_wide<unsigned __int128>(complete(70),true);
    check_wide<boost::multiprecision::uint256_t>(complete(240),true);
    check_wide<unsigned __int128>(split_graph(65,5),false);
    check_wide<boost::multiprecision::uint256_t>(split_graph(140,7,true),false);
    check_wide<boost::multiprecision::uint512_t>(split_graph(140,7,true),false);
    bool overflow=false;
    try{check_wide<uint64_t>(complete(70),true);}catch(const std::overflow_error&){overflow=true;}
    require(overflow,"missing count overflow rejection");
    std::cout<<"{\"passed\":true,\"graphs\":"<<graphs<<",\"checked_cells\":"<<cells
             <<",\"coverage_residual_checks\":"<<families<<",\"overflow_rejected\":true}\n";
}

static const char* names[]={"fixed","stream","active","plain","full","partial","replay-full","replay-partial"};
static constexpr int mode_count=sizeof(names)/sizeof(names[0]);

template<class T> struct Run {
    typename Kernel<T>::Data data;
    terminal::Metrics metrics;
    terminal::BuildWork construction;
    double build_ms=0,prepare_ms=0,choose_ms=0,solve_ms=0,cpu_ms=0,wall_ms=0;
    size_t index_bytes=0,paths=0,incidences=0,replay_bytes=0;
    uint64_t accepted=0,rejected=0;
};

template<class T> static Run<T> run(const Input& input,int maximum,int mode) {
    Run<T> result;const auto begin=Clock::now();const auto cpu=std::clock();
    if(mode<3) {
        auto start=Clock::now();Layout index(input.graph,maximum);result.build_ms=ms(start);
        start=Clock::now();index.prepare(input.graph.n);result.prepare_ms=ms(start);
        start=Clock::now();typename Kernel<T>::Combinations choose(input.d+1,maximum);result.choose_ms=ms(start);
        result.paths=index.paths.size();result.incidences=index.paths.vertices.size();result.index_bytes=index.bytes();
        start=Clock::now();
        if(mode==0)result.data=Kernel<T>{}.template fixed_sparse<true>(index,input.graph.n,choose,input.ordinary);
        else if(mode==1)result.data=std::move(Kernel<T>{}.template solve<true>(input.graph,index,choose,input.ordinary).common.data);
        else {
            auto replay=activeorder::Solver<T>{}.template solve<false>(input.graph,index,choose,input.ordinary,2);
            result.data=std::move(replay.common.data);result.accepted=replay.replay.accepted;result.rejected=replay.replay.rejected;
            result.replay_bytes=replay.replay.scratch_bytes;
        }
        result.solve_ms=ms(start);
    }else {
        auto start=Clock::now();terminal::Index index(maximum);terminal::build(input.graph,index,mode<6?mode-3:mode-5);result.build_ms=ms(start);
        start=Clock::now();index.prepare(input.graph.n);result.prepare_ms=ms(start);
        start=Clock::now();typename Kernel<T>::Combinations choose(input.d+1,maximum);result.choose_ms=ms(start);
        result.paths=index.rows.size();result.incidences=index.members.size();result.index_bytes=index.bytes();result.construction=index.work;
        start=Clock::now();auto output=mode<6?terminal::Solver<T>::solve(input.graph,index,choose,input.ordinary)
            :terminal::Solver<T>::template solve<false,true>(input.graph,index,choose,input.ordinary);result.solve_ms=ms(start);
        result.accepted=output.replay.accepted;result.rejected=output.replay.rejected;result.replay_bytes=output.replay.scratch_bytes;
        result.metrics=output.metrics;result.data=std::move(output.common.data);
    }
    result.wall_ms=ms(begin);result.cpu_ms=1000.0*(std::clock()-cpu)/CLOCKS_PER_SEC;return result;
}

template<class T> static void write(const Run<T>& result,int mode,int trial,int position) {
    const auto& work=result.data.work;const auto& metrics=result.metrics;const auto& build=result.construction;
    const double compute=result.build_ms+result.prepare_ms+result.choose_ms+result.solve_ms;
    std::cout<<"{\"mode\":\""<<names[mode]<<"\",\"trial\":"<<trial<<",\"position\":"<<position
        <<",\"build_ms\":"<<result.build_ms<<",\"prepare_ms\":"<<result.prepare_ms
        <<",\"choose_ms\":"<<result.choose_ms<<",\"solve_ms\":"<<result.solve_ms<<",\"compute_ms\":"<<compute
        <<",\"call_wall_ms\":"<<result.wall_ms<<",\"call_cpu_ms\":"<<result.cpu_ms
        <<",\"paths\":"<<result.paths<<",\"incidences\":"<<result.incidences<<",\"index_bytes\":"<<result.index_bytes
        <<",\"state_bytes\":"<<result.data.state_bytes<<",\"output_bytes\":"<<result.data.core.capacity()*sizeof(T)
        <<",\"source_reads\":"<<work.source_reads<<",\"target_reads\":"<<work.target_reads<<",\"count_reads\":"<<work.count_reads
        <<",\"events\":"<<work.events<<",\"heap_updates\":"<<work.updates<<",\"positive_writes\":"<<metrics.positive_writes
        <<",\"affected_rows\":"<<metrics.affected_rows<<",\"group_updates\":"<<metrics.group_updates
        <<",\"lazy_discards\":"<<metrics.lazy_discards<<",\"lazy_writes\":"<<metrics.lazy_writes
        <<",\"coefficient_calls\":"<<metrics.coefficient_calls<<",\"full_groups\":"<<build.full_groups
        <<",\"partial_groups\":"<<build.partial_groups<<",\"choice_members\":"<<build.choices
        <<",\"saved_members\":"<<build.saved_members<<",\"search_states\":"<<build.states
        <<",\"degree_tests\":"<<build.degree_tests<<",\"minimum_writes\":"<<build.minimum_writes
        <<",\"minimum_scratch_bytes\":"<<build.minimum_scratch_bytes<<",\"replay_bytes\":"<<result.replay_bytes
        <<",\"accepted\":"<<result.accepted<<",\"rejected\":"<<result.rejected<<'}';
}

template<class T> static void paired(const Input& input,int maximum,int rounds) {
    const auto expected=run<T>(input,maximum,0).data.core;
    for(int mode=0;mode<mode_count;++mode)require(run<T>(input,maximum,mode).data.core==expected,"paired warmup output mismatch");
    std::cout<<",\"core_hash\":\""<<hash_core(expected)<<"\",\"samples\":[";
    bool first=true;
    for(int trial=0;trial<rounds;++trial)for(int position=0;position<mode_count;++position) {
        const int mode=(trial+position)%mode_count;
        const auto result=run<T>(input,maximum,mode);
        require(result.data.core==expected,"paired timed full output mismatch");
        if(!first)std::cout<<',';first=false;write(result,mode,trial,position);
    }
    std::cout<<"],\"passed\":true}\n";
}

template<class T> static void single(const Input& input,int maximum,int mode) {
    const auto result=run<T>(input,maximum,mode);
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    size_t rss=usage.ru_maxrss;
#ifndef __APPLE__
    rss*=1024;
#endif
    std::cout<<",\"core_hash\":\""<<hash_core(result.data.core)<<"\",\"peak_rss_bytes\":"<<rss<<",\"sample\":";
    write(result,mode,0,0);std::cout<<",\"passed\":true}\n";
}

int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--selftest"){selftest();return 0;}
        require(argc>=4,"usage: terminal --paired graph rounds [max_s] | --single graph mode bits max_s");
        const std::string action=argv[1];const auto input=prepare(argv[2]);
        const bool is_single=action=="--single";
        require(is_single || action=="--paired","unknown action");
        int maximum=std::max(2,static_cast<int>(input.d)+1),mode=0,rounds=0;unsigned bits=0;
        if(is_single) {
            require(argc==6,"single needs mode, bits, max_s");
            while(mode<mode_count && names[mode]!=std::string(argv[3]))++mode;
            require(mode<mode_count,"unknown single mode");bits=std::stoul(argv[4]);maximum=std::stoi(argv[5]);
        }else {
            rounds=std::stoi(argv[3]);require(rounds>0,"positive trials required");if(argc>4)maximum=std::stoi(argv[4]);
            Layout index(input.graph,maximum);index.prepare(input.graph.n);bits=width(count_bound(input.graph,index,input.d));
        }
        require(maximum>=2 && maximum<=std::max(2,static_cast<int>(input.d)+1),"invalid maximum s");
        require(bits==64 || bits==128 || bits==256 || bits==512,"invalid count width");
        std::cout<<std::fixed<<std::setprecision(6)<<"{\"n\":"<<input.graph.n<<",\"m\":"<<input.graph.m
                 <<",\"maximum_s\":"<<maximum<<",\"ordinary_max\":"<<input.d<<",\"count_bits\":"<<bits
                 <<",\"load_ms\":"<<input.load_ms<<",\"order_ms\":"<<input.order_ms;
        auto dispatch=[&]<class T>(){if(is_single)single<T>(input,maximum,mode);else paired<T>(input,maximum,rounds);};
        if(bits==64)dispatch.template operator()<uint64_t>();
        else if(bits==128)dispatch.template operator()<unsigned __int128>();
        else if(bits==256)dispatch.template operator()<boost::multiprecision::uint256_t>();
        else dispatch.template operator()<boost::multiprecision::uint512_t>();
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
