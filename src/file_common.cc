#include <sys/stat.h>

int file_exists(const char * fn)
{
  struct stat buf;
  return (stat(fn,&buf) != -1);
}
