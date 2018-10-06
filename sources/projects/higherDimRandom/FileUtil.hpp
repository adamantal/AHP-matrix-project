#ifndef FILEUTIL_HPP
#define FILEUTIL_HPP

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <string>

std::set<std::string> scanFolder(std::string folder) {
    std::set<std::string> ret;
    DIR* dirp = opendir(folder.c_str());
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        ret.insert(dp->d_name);
    }
    (void)closedir(dirp);
    for (auto it = ret.begin(); it != ret.end(); ) {
        if (it->front() == '.') {
            ret.erase(it++);
        } else {
            ++it;
        }
    }
    return ret;
}

#endif
