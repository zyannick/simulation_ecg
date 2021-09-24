#pragma once


#include <iostream>
#include <cmath>
#include <experimental/filesystem>

using namespace std;




inline string separator()
{
#ifdef _WIN32
    return "\\";
#else
    return "/";
#endif
}

const string kPathSeparator =
        #ifdef _WIN32
        "\\";
#else
        "/";
#endif


std::vector<string> split(const string& str, const string& delim)
{
    std::vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

string concat_strings(std::vector<string>  vect_string, string delim)
{
    string new_string = "";
    for (size_t i =0; i < vect_string.size(); i++)
    {
        new_string += vect_string[i] + delim;
    }
    return new_string;
}


std::vector<string> get_subset(std::vector<string>  vect_string, size_t debut, size_t fin)
{
    std::vector<string> result;
    for(size_t i = debut; i < fin ; i++)
    {
        result.push_back(vect_string[i]);
    }
    return result;
}






bool create_directory(string dir_name, bool verbose = false)
{
    namespace fs = std::experimental::filesystem;
    string string_dir_name(dir_name);

    if(string_dir_name.size() < 1)
    {
        return false;
    }


    if(!fs::is_directory(string_dir_name) || !fs::exists(string_dir_name))
    {
        if(fs::create_directory(string_dir_name))
        {
            if(verbose)
                cout << "Directory " << string_dir_name << " created" << endl;
            return true;
        }
        else {
            if(verbose)
                cout << "Directory " << string_dir_name << " not created" << endl;
        }
    }
    else {
        if(verbose)
            cout <<  string_dir_name << " already exists or is not a directory" << endl;
    }

    return false;

}



bool create_directories(string dir_name, bool verbose = false)
{
    namespace fs = std::experimental::filesystem;
    string string_dir_name(dir_name);

    std::vector<string> list_dirs = split(string_dir_name, kPathSeparator);

    for(size_t i = 0; i < list_dirs.size(); i++)
    {
        string new_dir = concat_strings(get_subset(list_dirs, 0 , i+1), kPathSeparator);
        create_directory(new_dir, verbose);
    }


    return false;

}

void delete_before_create_directory(string dir_name)
{
    namespace fs = std::experimental::filesystem;
    string string_dir_name(dir_name);
    if(fs::is_directory(string_dir_name) && fs::exists(string_dir_name))
    {
        fs::remove_all(string_dir_name);
    }

    if(!fs::is_directory(string_dir_name) || !fs::exists(string_dir_name))
    {
        fs::create_directory(string_dir_name);
    }
}


string get_number_of_zero(int value, int number_of_zeros = 6)
{
    string nb_zeros = "";
    if(value > 0)
    {
        number_of_zeros = number_of_zeros - int(floor(log10(value))) ;
    }
    else
    {
        number_of_zeros = number_of_zeros - int(floor(log10(value+1))) ;
    }
    for(int i = 0; i < number_of_zeros; i++)
        nb_zeros = nb_zeros + "0";

    return nb_zeros;
}


