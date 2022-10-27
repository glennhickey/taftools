#include "uce.hpp"

#include <iostream>
#include <string>

using namespace std;


int main(int argc, char** argv) {

    string subcommand = argc > 1 ? argv[1] : "";
    if (subcommand != "uce") {
        cerr << "taftools: do something with taf files" << endl << endl
             << "usage: taftools <command> [options]" << endl << endl
             << "commands:" << endl
             << " -- uce     compute ultraconserved elements" << endl
             << endl;
        return 1;        
    }

    int ret = 1;
    if (subcommand == "uce") {
        ret = uce_main(argc - 1, argv + 1);
    }

    return ret;
}
