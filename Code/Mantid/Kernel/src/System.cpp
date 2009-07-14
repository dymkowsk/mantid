//-----------------------------
//Includes
//-----------------------------
#include "MantidKernel/System.h"
#include "Poco/Path.h"
#include <climits>
#include <cfloat>

// Need OS defined functions
#ifdef _WIN32
  #include <windows.h>
#elif defined __linux__
  #include <unistd.h>
  #include <fstream>
  #include <sstream>
  #include <algorithm>
  #include <iomanip>
  #include <iostream>
#elif defined __APPLE__
  #include <mach-o/dyld.h>
#endif

/**
 * Return what we consider to be an empty integer
 */
int Mantid::EMPTY_INT()
{
  return INT_MAX;
}

/**
 * Return what we consider to be an empty double
 */
double Mantid::EMPTY_DBL()
{
  return /*  -  */  DBL_MAX/2;
}

/**
 * Get the directory containing the program executable
 * @returns A string containing the path of the directory 
 * containing the executable, including a trailing slash
 */
std::string Mantid::Kernel::getDirectoryOfExecutable()
{
  return Poco::Path(getPathToExecutable()).parent().toString();
}

/**
  * Get the  full path the the program executable. This is not necessarily MantidPlot or a main
  * program since if we are running through Python, it is the Python executable that is considered
  * as the executing program
  * @returns A string containing the full path the the executable
  */
std::string Mantid::Kernel::getPathToExecutable()
{
  std::string execpath("");
  const size_t LEN(1024);
  char pBuf[LEN];
// The linux function returns an int, the Windows & Mac ones an unsigned type
#ifndef __linux__
  unsigned 
#endif
  int bytes(0);
  
#ifdef _WIN32
  bytes = GetModuleFileName(NULL, pBuf, LEN);
#elif defined __linux__
  char szTmp[32];
  sprintf(szTmp, "/proc/%d/exe", getpid());
  bytes = readlink(szTmp, pBuf, LEN);
#elif defined __APPLE__
  // Two calls to _NSGetExecutablePath required - first to get size of buffer
  _NSGetExecutablePath(pBuf,&bytes);
  const int success = _NSGetExecutablePath(pBuf,&bytes);
  if (success < 0) bytes = 1025;
#endif

  if( bytes > 0 && bytes < 1024 )
  {
    pBuf[bytes] = '\0';
    execpath = std::string(pBuf);
  }
  return execpath;
}

/**
 * Check if the path is on a network drive
 * @param path The path to be checked
 * @return True if the path is on a network drive.
 */
bool Mantid::Kernel::isNetworkDrive(const std::string & path)
{
#ifdef _WIN32
  // if path is relative get the full one
  char buff[MAX_PATH];
  GetFullPathName(path.c_str(),MAX_PATH,buff,NULL);
  std::string fullName(buff);
  size_t i = fullName.find(':');

  // if the full path doesn't contain a drive letter assume it's on the network
  if (i == std::string::npos) return true;

  fullName.erase(i+1);
  fullName += '\\';  // make sure the name has the trailing backslash
  UINT type = GetDriveType(fullName.c_str());
  return DRIVE_REMOTE == type;
#elif defined __linux__
  // This information is only present in the /proc/mounts file on linux. There are no drives on
  // linux only mount locations therefore the test will have to check the path against
  // entries in /proc/mounts to see if the filesystem type is NFS or SMB (any others ????)
  // Each line corresponds to a particular mounted location
  // 1st column - device name
  // 2nd column - mounted location
  // 3rd column - filesystem type commonly ext2, ext3 for hard drives and NFS or SMB for
  //              network locations

  std::ifstream mntfile("/proc/mounts");
  std::string txtread("");
  while( getline(mntfile, txtread) )
  {
    std::istringstream strm(txtread);
    std::string devname(""), mntpoint(""), fstype("");
    strm >> devname >> mntpoint >> fstype;
    if( !strm ) continue;
    // I can't be sure that the file system type is always lower case
    std::transform(fstype.begin(), fstype.end(), fstype.begin(), toupper);
    // Skip the current line if the file system isn't a network one
    if( fstype != "NFS" && fstype != "SMB" ) continue;
    // Now we have a line containing a network filesystem and just need to check if the path
    // supplied contains the mount location. There is a small complication in that the mount
    // points within the file have certain characters transformed into their octal 
    // representations, for example spaces->040.
    std::string::size_type idx = mntpoint.find("\\0");
    if( idx != std::string::npos ) 
    {
      std::string oct = mntpoint.substr(idx + 1, 3);
      strm.str(oct);
      int printch(-1);
      strm.setf( std::ios::oct, std::ios::basefield );  
      strm >> printch;
      if( printch != -1 )
      { 
        mntpoint = mntpoint.substr(0, idx) + static_cast<char>(printch) + mntpoint.substr(idx + 4);
      }
      // Search for this at the start of the path
      if( path.find(mntpoint) == 0 ) return true;
    }     
  }
  return false;
#else
  // Not yet implemented for the mac
  return false;
#endif
}
