#include <sys/stat.h>
#include <fstream>

using namespace std;

int file_exists(const char * fn)
{
  struct stat buf;
  return (stat(fn,&buf) != -1);
}

int is_gzip(const char * fn)
{
  ifstream in(fn);
  if(!in) return 0;
  int x;
  in.read(reinterpret_cast<char*>(&x),sizeof(int));
  auto first = x & 0xFF;
  auto second = ((x >> (8*1)) & 0xFF);
  return (first==0x1f) && (second==0x8b);
}
