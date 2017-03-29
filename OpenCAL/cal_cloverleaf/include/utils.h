#ifndef UTILS_H_
#define UTILS_H_
#include<string>



using std::string;


inline bool checkBoundaries(int x, int y, const int x_lb, const int x_up, const int y_lb, const int y_up){
    if(x < x_lb || x >= x_up || y < y_lb || y >= y_up)
        return false;
    return true;
}


FILE* open_FILE_orDie(const string& file_path, const string& modifiers){
    FILE * f;
    if ((f = fopen(file_path.c_str(),modifiers.c_str())) != NULL) {
        return f;
    }else{
        fprintf(stderr,"can't open file %s\nEXITING . . .\n",file_path.c_str());
        exit(-1);
    }
}

#endif //UTILS_H_
