//
// Daq common tools
//
//

#ifndef __DQERROR__H
#define __DQERROR__H

#include <cstdint>
#include <bit>
#include <TSysEvtHandler.h>
#include <TString.h>
#include <Rtypes.h>

// DQ_FileTest
#include <dirent.h>
#include <fnmatch.h>
#include <Riostream.h>
#include <TString.h>
#include <TObjArray.h>

//#define uint32_t Int_t
//#define uint16_t Short_t

// --- common parameters

#define DQ_LASTRUN_FILE "../cfg/last.run"

// --- verbosity macro

#define DQ_VERBOSE 3 // 0= information, 1= +errror, 2= +warning, 3= +message, 4= +debug message, 5= +dump 

#define DQ_STR "%s[%d] "
#define DQ_ARG __FUNCTION__, __LINE__
#define DQ_PRINT(...) fprintf(stderr,__VA_ARGS__)

#if (DQ_VERBOSE>4)
#define DQ_DUM(_fmt, ...) DQ_PRINT("DUM:" DQ_STR _fmt, DQ_ARG, ##__VA_ARGS__)
#else
#define DQ_DUM(...)
#endif
#if (DQ_VERBOSE>3)
#define DQ_DBG(_fmt, ...) DQ_PRINT("DBG:" DQ_STR _fmt, DQ_ARG, ##__VA_ARGS__)
#else
#define DQ_DBG(...)
#endif
#if (DQ_VERBOSE>2)
#define DQ_MSG(_fmt, ...) DQ_PRINT("MSG:" DQ_STR _fmt, DQ_ARG, ##__VA_ARGS__)
#else
#define DQ_MSG(...)
#endif
#if (DQ_VERBOSE>1)
#define DQ_WRN(_fmt, ...) DQ_PRINT("WRN:" DQ_STR _fmt, DQ_ARG, ##__VA_ARGS__)
#else
#define DQ_WRN(...)
#endif
#if (DQ_VERBOSE>0)
#define DQ_ERR(_fmt, ...) DQ_PRINT("ERR:" DQ_STR _fmt, DQ_ARG, ##__VA_ARGS__)
#else
#define DQ_ERR(...)
#endif
#define DQ_INF(_fmt, ...) DQ_PRINT("INF:" DQ_STR _fmt, DQ_ARG, ##__VA_ARGS__)

#define GEM_MAXERROR 64

#define GEM_ERR_APVNUMBER 0
#define GEM_ERR_TRAILER 1
#define GEM_ERR_HEADER 2
#define GEM_ERR_HADDR 3
#define GEM_ERR_HERR 4
#define GEM_ERR_SAMPLENUMBER 5
#define GEM_ERR_WCOUNT 6
#define GEM_ERR_DATA 7
#define GEM_ERR_EBLOCK 8
#define GEM_ERR_TAG 9

class DQerror {

 public:

  DQerror();

  ~DQerror();

  void inc(Int_t imask, Int_t step=1) {
    fErrorCount[imask] += step;
    fErrorMask |= (1<<imask);
  };

  ULong64_t getMask() { return fErrorMask; };
  void resetMask(ULong64_t v=0) { fErrorMask = v; };

  void print(Float_t norm=1.);

 private:

  Int_t fErrorCount[GEM_MAXERROR];
  ULong64_t fErrorMask;
  TString fErrorName[GEM_MAXERROR];

  Int_t fErrorMax;
};


class DQSigHandler : public TSignalHandler {
 public:
 DQSigHandler() : TSignalHandler(kSigInterrupt, kFALSE) { 
    fctrlc_pressed=0;
  }

  virtual Bool_t Notify()
  {
    fctrlc_pressed++;
    DQ_MSG("CTRL-C pressed, handled by DQSigHandler %d\n",fctrlc_pressed);
    return kTRUE;

  }
  
  Int_t Ctrlc_IsPressed() { return fctrlc_pressed; }

 private:
  Int_t fctrlc_pressed;

  ClassDef (DQSigHandler, 0) 

};

/**
 * File Utility Class
 * 
 */

class DQFile {

 private:

  Int_t fLastRunMode; // =1 if use lastrun file content

 public:

  DQFile() {

    fLastRunMode=0;
  };

  ~DQFile() { };

    /**     
     Try to read file with last run number
    */
    
  int ReadLastRunNumber(const char *runfile) {
    FILE *ff;
    int rr;
    ff = fopen(runfile,"r");
    if (ff) {
      int no = fscanf(ff,"%d",&rr);
      fclose(ff);
    } else {
      rr = -1;
    }
    return rr;
  };

  /**
     Try to match file pattern in directory
     iDir: search dir
     iFile: file name pattern searched in iDir
     
     return full path on success, empty string if nothing found
  */

  TString FileMatchInDir(TString iDir, TString iFile) {

    TString gfile;

    DIR *dir = opendir(iDir.Data());
    if (dir) {
      struct dirent *ent;
      while ((ent = readdir(dir)) != NULL) {
	if (fnmatch(iFile.Data(), (ent->d_name),  FNM_PATHNAME) == 0) {
	  gfile = iDir + (ent->d_name); 
	  DQ_MSG("Found file %s\n", gfile.Data());
	  return gfile;
	}
      }
    }
    
    DQ_WRN("No file %s found in dir %s\n", iFile.Data(), iDir.Data());
    return "";

  };

  /** 
    Try to discover file input looking at iPath + iFile
    or tring to read last run number from iFile and than
    looking at iPath + iFile
    
    iPath: search directory path
    iFile: either input file (search first) or full path of last.run file (second choice)
    irunPattern: the generic pattern of the file name with run number
    Return: full file path
    : empty string if nothing found
    
  */
  
  TString FileDiscovery(TString iPath, TString iFile, TString irunPattern="*%d.dat") {
    
    Int_t lrun=-2;
    // try file : iPath + iFile
    TString gfile;

    if (fLastRunMode==0) {

      DQ_MSG("Try input data from %s / %s\n",iPath.Data(),iFile.Data()); 

      gfile = FileMatchInDir(iPath, iFile);
    
      if (gfile.Length()>0) { 
	return gfile;
      }
    
      //  DQ_MSG("File %s/%s not found, try last run data\n",iPath.Data(),iFile.Data());
    
    }

    DQ_MSG("Try to read last run number from %s \n",iFile.Data()); 

    lrun = ReadLastRunNumber(iFile.Data());
    if (lrun == -1) {
      DQ_ERR("Cannot open last run file %s\n",iFile.Data());
      return "";
    }

    fLastRunMode=1;
    
    // try last run file
    gfile = FileMatchInDir(iPath, Form(irunPattern.Data(),lrun));
    
    if (gfile.Length()>0) {
      DQ_MSG("***** Got file from last run number %d *****\n",lrun);
      return gfile;
    }
    
    DQ_WRN("Not able to open input file from last run number %d, maybe too early\n",lrun);
    
    return "__WAIT__";

  };

  Int_t LastRunMode() {
    return fLastRunMode;
  };
  
};

// Other utility methods

// strip string and append string

TString DQ_SandA(const TString &s1,  // original string
		 const TString &sa,  // appended string
		 const char c='.')   // strip from last 'c' character (included) 

{
  Ssiz_t iop = s1.Last(c);
  if (iop == kNPOS) { iop = s1.Length()-1; }
  TString so = s1(0,iop) + sa;
  return so;
};


/*
 * Triangular Array Linearization
 *
 * Array of n*n elements -> triangular (including diagonal) = n*(n+1)/2
 *
 *                   i  
 *                   V
 *    0  1  3  6 10  .  . 
 *       2  4  7 11  .  .
 *  j>      5  8 12  k  .
 *             9 13  .  .
 *               14  .  .
 *                   .  .
 *                      .
 *
 *   k = i*(i+1)/2 + j 
 */

template <class T> class DQ_TriArr {
 private:
  T* arr;
  Int_t fN;
 public:
  DQ_TriArr(Int_t n) {
    fN = n;
    arr = new T[n*(n+1)/2];
  };
  ~DQ_TriArr() {
    delete[] arr;
  };
  T Get(Int_t i, Int_t j) {
    return (T) arr[i*(i+1)/2+j];
  };
  void Set(Int_t i, Int_t j, T v) {
    arr[i*(i+1)/2+j]=v;
  };
  Int_t Size() { return fN; };
  Int_t Index(Int_t i, Int_t j) { return i*(i+1)/2+j; };

  T* Data() { return (T*) arr; };
};

#endif
