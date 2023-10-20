/*
 * Source Combiner
 * @Author: Rami Zouari
 * This program takes a source file and combines all the included files into one file
 * This is useful for submitting a single file to online judges.
 *
 * @Usage:
 * source_combiner -i <input file> -o <output file> -p <paths to search for files> [-v -e]
 * If no input file is specified, the program reads from stdin
 * If no output file is specified, the program writes to stdout
 *
 * @Options:
 * -i <input file> : The input file to read from
 * -o <output file> : The output file to write to
 * -p <paths to search for files> : The paths to search for files
 * -v : Verbose output
 * -e : Force equivalence of the paths. This will guarantee that each file is included only once.
 *
 * @Notes:
 * The program will replace the first occurrence of #include "file" with the contents of file. Note that this is not a full C++ preprocessor, and
 * it may output a wrong result if the input file contains #include directives in comments or strings.
 *
 * @Warning
 * The author is not responsible for any damage caused by this program. Use at your own risk.
 * */

#include <iostream>
#include <filesystem>
#include <utility>
#include <vector>
#include <fstream>
#include "boost/program_options.hpp"
#include <regex>
#include <unordered_set>
#include <optional>

class Paths
{
    std::vector<std::filesystem::path> paths;
public:
    explicit Paths(std::vector<std::string> paths) : paths(paths.begin(),paths.end())
    {
        paths.emplace_back(std::filesystem::current_path());
    }

    explicit Paths(const std::string& representation) : Paths()
    {
        std::stringstream ss(representation);
        std::string path;
        while(std::getline(ss,path,':'))
            paths.emplace_back(path);
    }

    Paths():paths{std::filesystem::current_path()}
    {
    }

    std::filesystem::path find(const std::filesystem::path &file,const std::filesystem::path &base)
    {
        auto f = try_find(file,base);
        if(f) return *f;
        else throw std::runtime_error("File " + file.string() + " Not found");
    }

    std::filesystem::path find(const std::filesystem::path &file)
    {
        return find(file,"");
    }

    std::optional<std::filesystem::path> try_find(const std::filesystem::path &file,const std::filesystem::path &base)
    {
        auto candidate=base/file;
        if(std::filesystem::exists(candidate))
            return candidate;
        for(auto& p:paths)
        {
            candidate = p / file;
            if (std::filesystem::exists(candidate))
                return candidate;

        }
        return std::nullopt;
    }

    std::optional<std::filesystem::path> try_find(const std::filesystem::path &file)
    {
        return try_find(file,"");
    }

    Paths& operator+=(const Paths &other)
    {
        paths.insert(paths.end(),other.paths.begin(),other.paths.end());
        return *this;
    }

    Paths operator+(const Paths &other) const
    {
        Paths p(*this);
        p+=other;
        return p;
    }

    Paths& operator+=(const std::filesystem::path &other)
    {
        paths.push_back(other);
        return *this;
    }

    Paths operator+(const std::filesystem::path &other) const
    {
        Paths p(*this);
        p+=other;
        return p;
    }

    Paths& operator+=(const std::string &other)
    {
        return *this+=Paths(other);
    }
};

class strong_path : public std::filesystem::path
{
    using std::filesystem::path::path;
    inline static Paths *paths;
    std::filesystem::path base;
public:
    static void register_paths(Paths &_path)
    {
        paths = &_path;
    }
    strong_path(const std::filesystem::path &p): std::filesystem::path(p){}
    strong_path(const std::filesystem::path &p, std::filesystem::path base): std::filesystem::path(p), base(std::move(base)){}
    strong_path(std::filesystem::path && p): std::filesystem::path(std::move(p)){}
    bool operator==(const strong_path &other) const
    {
        using parent_t=std::filesystem::path;
        auto x=paths->try_find(*this,base);
        auto y=paths->try_find(other,other.base);
        if(x && y)
            return std::filesystem::equivalent(*x,*y);
        return x==y;
    }

    [[nodiscard]] std::filesystem::path weak() const
    {
        return *this;
    }

};

class UnorderedUnionFind
{
    using WeakPathType=std::filesystem::path;
    using StrongPathType=strong_path;
    struct ParentData
    {
        const WeakPathType* parent;
        int rank;
    };
    std::unordered_map<WeakPathType,ParentData> data{};
    std::unordered_set<WeakPathType> representatives;
public:
    void connect(const StrongPathType& A, const StrongPathType&B)
    {
        if(equivalent(A,B))
            return;
        auto &[C,NodeData1]=get(representative(A));
        auto &[D,NodeData2]=get(representative(B));
        auto &[p1,r1]=NodeData1;
        auto &[p2,r2]=NodeData2;
        if(r1<r2) {
            p1 = p2;
            representatives.erase(C);
        }
        else if(r1>r2) {
            p2 = p1;
            representatives.erase(D);
        }
        else
        {
            p1=p2;
            r2++;
            representatives.erase(C);
        }
    }
    std::pair<const WeakPathType,ParentData>& get(const StrongPathType &A)
    {
        if(!data.contains(A))
        {
            auto [it,_]=data.emplace(A,ParentData{nullptr,0});
            auto &[B,NodeData]=*it;
            NodeData.parent=&B;
            bool has_representative=false;
            for(const auto &R: representatives) if (A == strong_path{R})
            {
                NodeData.parent = &representative(R);
                has_representative=true;
            }
            if(!has_representative)
                representatives.insert(A);
            return *it;
        }
        return *data.find(A);
    }
    const WeakPathType& representative(const StrongPathType&A)
    {
        auto &[B,NodeData]=get(A);
        auto &C=NodeData.parent;
        if(&B==C)
            return B;
        else
        {
            NodeData.parent = &representative(*NodeData.parent);
            return *NodeData.parent;
        }
    }

    bool equivalent(const StrongPathType&A, const StrongPathType&B)
    {
        return &representative(get(A).first)== &representative(get(B).first);
    }

};


template<>
struct std::hash<strong_path>
{
    std::size_t operator()(const strong_path &p) const
    {
        return std::filesystem::hash_value(p);
    }
};

class FileSet
{
    UnorderedUnionFind UF;
public:
    using equivalent_file_set_t = std::unordered_set<std::filesystem::path>;
    auto insert(const strong_path & path)
    {
        return files.insert(UF.get(path).first);
    }

    auto find(const strong_path & path)
    {
        return files.find(UF.representative(path));
    }

    auto begin() const
    {
        return files.begin();
    }

    auto end() const
    {
        return files.end();
    }

private:
    equivalent_file_set_t files;
};

struct Trie
{
    struct Node
    {
        std::unordered_map<char,Node> children;
        bool is_end=false;
    };
    Node root;
    void insert(const std::string &s)
    {
        auto node=&root;
        for(auto c:s)
        {
            auto it=node->children.find(c);
            if(it==node->children.end())
                it=node->children.emplace(c,Node{}).first;
            node=&it->second;
        }
        node->is_end=true;
    }

    bool contains(const std::string &s)
    {
        auto node=&root;
        for(auto c:s)
        {
            auto it=node->children.find(c);
            if(it==node->children.end())
                return false;
            node=&it->second;
        }
        return node->is_end;
    }

    bool contains(const char *w)
    {
        auto node=&root;
        for(auto c=w;*c;c++)
        {
            auto it=node->children.find(*c);
            if(it==node->children.end())
                return false;
            node=&it->second;
        }
        return node->is_end;
    }
};



std::vector<std::string> standard_cxx_headers {
        "iostream","concepts","chrono","compare","complex","exception",
        "format","functional","iterator","limits","list","map","memory","new","numeric","numbers","optional","queue","random",
        "ratio","regex","set","sstream","stdexcept","string","string_view","tuple","type_traits","typeindex", "typeinfo",
        "unordered_map","unordered_set","utility","valarray","vector", "filesystem","csetjmp","csignal","cstdarg","cstddef",
        "cstdio","cstdlib", "ctime","expected","cuchar","cwchar","cwctype","iosfwd","ios","istream","ostream","fstream",
        "strstream","iomanip", "initializer_list","locale","streambuf","cinttypes","cstdint","cfloat","climits","cstdbool",
        "cctype","cstring", "optional","variant","any","bitset","deque","forward_list","list","map","queue","set","stack",
        "source_location","version","memory_resource","scoped_allocator","execution","stop_token","semaphore","memory",
        "stdfloat","cassert","cerrno","stacktrace","system_error","charconv","array","flat_map","flat_set","span","mdspan",
        "generator","ranges","syncstream","future","mutex","shared_mutex","thread","latch","barrier","atomic",
        "condition_variable","bit","bitops","cfenv","clocale", "codecvt", "text_encoding","print","spanstream",
        "hazard_pointer","rcu","ciso646","cstdalign","cmath","algorithm"
};

std::vector<std::string> standard_c_headers {
        "assert.h", "ctype.h", "errno.h", "fenv.h", "float.h", "limits.h", "locale.h", "math.h", "setjmp.h", "signal.h",
        "stdarg.h", "stdbool.h", "stddef.h", "stdint.h", "stdio.h", "stdlib.h", "string.h", "tgmath.h", "threads.h",
        "inttypes.h","time.h","uchar.h","wchar.h","wctype.h","stdatomic.h","iso646.h","stdalign.h"
};

class Combiner
{
    Paths  &paths;
    Trie exclusion_trie;
    std::pair<bool,std::string> match_include_directive(const std::string &d)
    {
        std::regex r(R"~~(\s*#include\s*"(.+)".*)~~");
        std::smatch m;
        if(std::regex_match(d,m,r))
            return {true,m[1]};
        else if(accept_angle)
        {
            r=R"~~(\s*#include\s*<(.+)>.*)~~";
            if(std::regex_match(d,m,r) && ! exclusion_trie.contains(m[1]))
                return {true,m[1]};
        }
        return {false,""};
    };

    bool is_white_listed(const std::string &X)
    {
        return !white_list_regex.has_value() || std::regex_match(X,*white_list_regex);
    }

    using equivalent_file_set_t = FileSet;
    using equal_representation_set_t = std::set<std::filesystem::path>;

    template<typename FileSet>
    void combine(std::istream &in, std::ostream &out,FileSet &visitedFiles,const std::filesystem::path &base = std::filesystem::current_path())
    {
        std::string line;
        while(std::getline(in,line))
        {
            auto [is_include,file] = match_include_directive(line);
            if(is_include && is_white_listed(file))
            {
                std::filesystem::path filePath;
                if(ignore_not_exists)
                {
                    auto filePathOpt=paths.try_find(file,base);
                    if(!filePathOpt)
                    {
                        out << line << '\n';
                        continue;
                    }
                    else filePath=*filePathOpt;
                }
                else filePath=paths.find(file,base);
                auto weakPath=strong_path(filePath, base);
                auto it=visitedFiles.find(weakPath);
                if(it!=visitedFiles.end())
                    continue;
                visitedFiles.insert(weakPath);
                std::ifstream includedFile(filePath);
                combine(includedFile,out,visitedFiles,filePath.parent_path());
            }
            else
                out << line << '\n';
        }
    }
    bool accept_angle;
    bool ignore_not_exists;
    std::optional<std::regex> white_list_regex;
public:
    enum Flags : unsigned int
    {
        None = 0,
        AcceptAngle = 0b1,
        IgnoreNotExists = 0b10,
    };

    enum PathEquivalenceCriteria
    {
        SameRepresentation,
        SameFile,
    };

    Combiner(Paths &paths,Flags flags= None, const std::vector<std::string>& additional_excluded_files={}): paths(paths), accept_angle(flags & AcceptAngle),
        ignore_not_exists(flags & IgnoreNotExists)
    {
        for(auto &s:standard_c_headers)
            exclusion_trie.insert(s);
        for(auto &s:standard_cxx_headers)
            exclusion_trie.insert(s);
        for(auto &s:additional_excluded_files)
            exclusion_trie.insert(s);
    }

    void set_white_list_regex(const std::string& regex)
    {
        if(regex.empty()) white_list_regex=std::nullopt;
        else white_list_regex.emplace(regex);
    }

    void combine(std::istream &in, std::ostream &out, PathEquivalenceCriteria criteria=SameRepresentation)
    {
        if(criteria==SameRepresentation)
        {
            equal_representation_set_t visitedFiles;
            combine(in,out,visitedFiles);
        }
        else
        {
            equivalent_file_set_t visitedFiles;
            combine(in,out,visitedFiles);
        }

    }

};

Combiner::Flags operator|(Combiner::Flags a,Combiner::Flags b)
{
    return static_cast<Combiner::Flags>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
}

Combiner::Flags operator&(Combiner::Flags a,Combiner::Flags b)
{
    return static_cast<Combiner::Flags>(static_cast<unsigned int>(a) & static_cast<unsigned int>(b));
}

Combiner::Flags& operator|=(Combiner::Flags &a,Combiner::Flags b)
{
    a=static_cast<Combiner::Flags>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
    return a;
}

Combiner::Flags& operator&=(Combiner::Flags &a,Combiner::Flags b)
{
    a=static_cast<Combiner::Flags>(static_cast<unsigned int>(a) & static_cast<unsigned int>(b));
    return a;
}

int main(int argc,char **argv)
{

    namespace po = boost::program_options;
    po::options_description desc("Command Line options");
    unsigned int cpus=1;
    desc.add_options()
            ("help,h", "produce help message")
            ("output,o", po::value<std::filesystem::path>(), "Output file for the results")
            ("input,i", po::value<std::filesystem::path>(), "Input source file")
            ("verbose,v", po::bool_switch(), "Verbose output")
            ("paths,p", po::value<std::string>()->default_value(""), "Path to search for files")
            ("use-cpath,c", po::bool_switch(), "Use CPATH environment variable")
            ("ignore",po::value<std::string>()->default_value(""), "Ignore files matching this regex")
            ("brackets,b",po::bool_switch(),"Accept bracket includes. This is not recommended")
            ("ignore-not-exists",po::bool_switch(),"Ignore files that were not found")
            ("white-list,w",po::value<std::string>()->default_value(""),"White list regex")
            ("equivalent,e",po::bool_switch(), "File equivalence criteria");
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    std::ostream* outputStream;
    std::unique_ptr<std::ofstream> outFileStream;
    std::unique_ptr<std::ifstream> inFileStream;
    std::istream *inputStream;
    Paths paths(vm["paths"].as<std::string>());
    if(vm["use-cpath"].as<bool>())
    {
        auto cpath=std::getenv("CPATH");
        if(cpath)
            paths+=Paths(cpath);
    }
    strong_path::register_paths(paths);
    bool verbose=vm["verbose"].as<bool>();
    std::string begin_comment,end_comment;
    if(!vm.count("output"))
    {
        begin_comment = "/*\n";
        end_comment = "*/\n";
    }
    std::cout << begin_comment;
    if(!vm.count("input"))
        inputStream = &std::cin;
    else
    {
        auto file=vm["input"].as<std::filesystem::path>();
        inFileStream = std::make_unique<std::ifstream>(file);
        inputStream = inFileStream.get();
        if(verbose) std::cout << "Input file: " << file << std::endl;

    }

    if (!vm.count("output"))
        outputStream = &std::cout;
    else
    {
        auto file=vm["output"].as<std::filesystem::path>();
        outFileStream = std::make_unique<std::ofstream>(vm["output"].as<std::filesystem::path>());
        outputStream = outFileStream.get();
        if(verbose) std::cout << "Output file: " << file << std::endl;
    }
    std::vector<std::string> ignored_files;
    std::cout << end_comment;
    {
        std::stringstream ss(vm["ignore"].as<std::string>());
        std::string file;
        while(std::getline(ss,file,':'))
            ignored_files.emplace_back(file);
    }
    Combiner::Flags combiner_flags=Combiner::Flags::None;
    if(vm["brackets"].as<bool>())
        combiner_flags|=Combiner::Flags::AcceptAngle;
    if(vm["ignore-not-exists"].as<bool>())
        combiner_flags|=Combiner::Flags::IgnoreNotExists;
    Combiner combiner(paths,combiner_flags,ignored_files);
    combiner.set_white_list_regex(vm["white-list"].as<std::string>());
    std::set<std::filesystem::path> visitedFiles;
    if(vm["equivalent"].as<bool>())
        combiner.combine(*inputStream,*outputStream,Combiner::SameFile);
    else
        combiner.combine(*inputStream,*outputStream,Combiner::SameRepresentation);
}