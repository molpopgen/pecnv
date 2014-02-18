#!sh

find . -perm -111 -type f -maxdepth 1 -exec cp "{}" ~/bin \;