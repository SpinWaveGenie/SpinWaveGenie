/*
Copyright (C) 2011,2012 Remik Ziemlinski. See MIT-LICENSE.

CHANGELOG

v0.0.0 20110502 rsz Created.
V1.0.0 20110522 rsz Extended to show eta with growing bar.
v2.0.0 20110525 rsz Added time elapsed.
v2.0.0 20110525 rsz Fix minutes value.
v2.0.1 20111006 rsz Added default constructor value.
v2.0.2 20111011 rsz Fix minutes value (for real this time). Show decimal rate if less than one.
v2.0.3 20111109 rsz Fix overwrite previous line if longer than new line.
v2.1.0 20111130 rsz Templatized counter to make increments possible with reals. Replaced cout with faster printf.
*/

#ifndef EZ_RATEPROGRESSBAR_H
#define EZ_RATEPROGRESSBAR_H

#include <stdint.h>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <boost/format.hpp>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

namespace ez
{
/*
Display progress on two lines with rate update.

Done |  Elapsed | Remaining | Processed | Unprocessed | Rate
 90% | 23:45:56 |  23:45:56 |      1 MB |    1,299 MB | 124 MB/s
*/
template <typename T> class ezRateProgressBar
{
public:
  ezRateProgressBar(T _n = 0) : n(_n), cur(0), done(0), pct(0), width(0), unitsWidth(0) {}
  void reset()
  {
    pct = 0;
    cur = 0;
    done = 0;
  }
  void start()
  {
#ifdef WIN32
    assert(QueryPerformanceFrequency(&g_llFrequency) != 0);
#endif
    startTime = osQueryPerfomance();
    prevTime = startTime;
    printHeader();
  }

  void commaNumber(T v, std::string &str)
  {
    std::ostringstream stm;
    if (v >= 1)
    {
      stm << (long long)v;
      str = stm.str();
      for (int64_t i = str.size() - 3; i > 0; i -= 3)
        str.insert(i, 1, ',');
    }
    else
    {
      stm << std::fixed << std::setprecision(3) << v;
      str = stm.str();
    }
  }

  // http://stackoverflow.com/questions/3283804/c-get-milliseconds-since-some-date
  long long osQueryPerfomance()
  {
#ifdef WIN32
    LARGE_INTEGER llPerf = {0};
    QueryPerformanceCounter(&llPerf);
    return llPerf.QuadPart * 1000ll / (g_llFrequency.QuadPart / 1000ll);
#else
    struct timeval stTimeVal;
    gettimeofday(&stTimeVal, nullptr);
    return stTimeVal.tv_sec * 1000000ll + stTimeVal.tv_usec;
#endif
  }

  void printHeader()
  {
    // Done |  Elapsed | Remaining | Processed |  Unprocessed | Rate
    std::string out;
    out.reserve(80);
    out.assign("Done |  Elapsed | Remaining | ");
    // Compute how many max characters 'n' would take when formatted.
    std::string str;
    commaNumber(n, str);
    // Width of "123,456".
    std::size_t nWidth = str.size();
    unitsWidth = strlen(units);
    // Width of "123,456 MB".
    std::size_t nUnitsWidth = nWidth + unitsWidth + 1;
    // "Processed" is 9 chars.
    if (nUnitsWidth > 9)
      procColWidth = nUnitsWidth;
    else
      procColWidth = 9;

    if (nUnitsWidth > 11)
      unprocColWidth = nUnitsWidth;
    else
      unprocColWidth = 11;

    out.append(procColWidth - 9, ' ');
    out.append("Processed | ");
    out.append(unprocColWidth - 11, ' ');
    out.append("Unprocessed | Rate\n");

    printf("%s", out.c_str());
  }

  void secondsToString(unsigned int t, std::string &out)
  {
    unsigned int days = t / 86400;
    unsigned int sec = t - days * 86400;
    unsigned int hours = sec / 3600 + days * 24;
    sec = t - hours * 3600;
    sec = (sec > 0 ? sec : 0);
    unsigned int mins = sec / 60;
    sec -= mins * 60;
    std::stringstream tmp;
    out.clear();
    tmp << boost::format("%02u:%02u:%02u") % hours % mins % sec;
    out.append(tmp.str());
  }

  void operator++() { update(cur++); };

  T operator+=(const T delta)
  {
    cur += delta;
    update(cur);
    return cur;
  };

  void update(T newvalue)
  {
    // Nothing to update if already maxed out.
    if (done)
      return;
    endTime = osQueryPerfomance();

    // Abort if at least 1 second didn't elapse, unless newvalue will get us to 100%.
    if (((endTime - prevTime) / 1000000.0 < 1.0) && (newvalue < n))
      return;
    prevTime = endTime;
    double dt = (endTime - startTime) / 1000000.0;
    // if (dt < 1) return; // Was meant to avoid division by zero when time was in whole numbers.
    cur = newvalue;
    std::stringstream pctstr;
    double Pct = ((double)cur) / n;
    pctstr << boost::format("%3d%%") % static_cast<int>(100 * Pct);
    std::string out;
    out.reserve(80);
    out.append(pctstr.str());
    // Seconds.
    std::string tstr;
    out.append(" | ");
    secondsToString((unsigned int)dt, tstr);
    out.append(tstr);
    int64_t pad, newwidth;
    if (Pct >= 1.0)
    {
      // Print overall time and newline.
      out.append(" |  00:00:00 | ");
      commaNumber(n, tstr);
      pad = procColWidth - tstr.size() - unitsWidth - 1;
      if (pad > 0)
        out.append(pad, ' ');
      out.append(tstr);
      out.append(" ");
      out.append(units);
      out.append(" | ");
      pad = unprocColWidth - 1 - unitsWidth - 1;
      if (pad > 0)
        out.append(pad, ' ');
      out.append("0 ");
      out.append(units);
      out.append(" | ");
      commaNumber((unsigned int)(n / dt), tstr);
      out.append(tstr);
      out.append(" ");
      out.append(units);
      out.append("/s");

      newwidth = out.size();
      if (newwidth < width)
        out.append(width - newwidth, ' ');

      // Store length of this string so we know how much to blank out later.
      width = newwidth;
      out.append("\n");
      printf("%s", out.c_str());
      fflush(stdout);
      done = 1;
    }
    else
    {
      double eta = 0.;
      if (Pct > 0.0)
      {
        eta = dt * (1.0 - Pct) / Pct;
      }
      out.append(" |  ");
      secondsToString((unsigned int)eta, tstr);
      out.append(tstr);
      out.append(" | ");
      commaNumber(cur, tstr);
      pad = procColWidth - tstr.size() - unitsWidth - 1;
      if (pad > 0)
        out.append(pad, ' ');
      out.append(tstr);
      out.append(" ");
      out.append(units);
      out.append(" | ");
      commaNumber(n - cur, tstr);
      pad = unprocColWidth - tstr.size() - unitsWidth - 1;
      if (pad > 0)
        out.append(pad, ' ');
      out.append(tstr);
      out.append(" ");
      out.append(units);
      out.append(" | ");
      eta = cur / dt;
      if (eta > 1.0)
        commaNumber((unsigned int)eta, tstr);
      else
      {
        std::ostringstream stm;
        stm << eta;
        tstr = stm.str();
      }
      out.append(tstr);
      out.append(" ");
      out.append(units);
      out.append("/s");
      // Pad end with spaces to overwrite previous string that may have been longer.
      newwidth = out.size();
      if (newwidth < width)
        out.append(width - newwidth, ' ');

      width = newwidth;
      out.append("\r");
      printf("%s", out.c_str());
      fflush(stdout);
    }
  }

  T n;
  T cur;
  char done;
  unsigned short pct; // Stored as 0-1000, so 2.5% is encoded as 25.
  int64_t width;      // Length of previously printed line we need to overwrite.
  std::size_t unitsWidth;
  std::size_t procColWidth;
  std::size_t unprocColWidth;
  const char *units; // Unit string to show for processed data.
  long long startTime, prevTime, endTime;
#ifdef WIN32
  LARGE_INTEGER g_llFrequency;
#endif
};
}
#endif // EZ_RATEPROGRESSBAR_H
