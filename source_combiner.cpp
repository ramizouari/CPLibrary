/*
 * Source Combiner
 * @Author: Rami Zouari
 * This program takes a source file and combines all the included files into one file
 * This is useful for submitting a single file to online judges.
 *
 * @Usage:
 * source_combiner -i <input file> -o <output file> -p <paths to search for files>
 * If no input file is specified, the program reads from stdin
 * If no output file is specified, the program writes to stdout
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
#include <vector>
#include <fstream>
#include "boost/program_options.hpp"
#include <regex>

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
        auto candidate=base/file;
        if(std::filesystem::exists(candidate))
            return candidate;
        for(auto& p:paths)
        {
            candidate = p / file;
            if (std::filesystem::exists(candidate))
                return candidate;

        }
        throw std::runtime_error("File " + file.string() + " Not found");
    }

    std::filesystem::path find(const std::filesystem::path &file)
    {
        return find(file,"");
    }
};

class Combiner
{
    Paths  &paths;


    std::pair<bool,std::string> match_include_directive(const std::string &d)
    {
        std::regex r(R"~~(\s*#include\s*"(.+)".*)~~");
        std::smatch m;
        if(std::regex_match(d,m,r))
            return {true,m[1]};
        return {false,""};
    };

public:
    Combiner(Paths &paths):paths(paths)
    {
    }


    void combine(std::istream &in, std::ostream &out,std::set<std::filesystem::path> &visitedFiles,std::filesystem::path base = {})
    {
        std::string line;
        while(std::getline(in,line))
        {
            auto [is_include,file] = match_include_directive(line);
            if(is_include)
            {
                if(visitedFiles.find(file)!=visitedFiles.end())
                    continue;
                visitedFiles.insert(file);
                auto filePath=paths.find(file,base);
                std::ifstream includedFile(filePath);
                combine(includedFile,out,visitedFiles,filePath.parent_path());
            }
            else
                out << line << '\n';
        }
    }

};

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
            ("paths,p", po::value<std::string>(), "Path to search for files");
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
    std::cout << end_comment;
    Combiner combiner(paths);
    std::set<std::filesystem::path> visitedFiles;
    combiner.combine(*inputStream,*outputStream,visitedFiles);
}