#ifndef _ProgressBar_Class
#define _ProgressBar_Class
#include <iostream>
#include <string>
#include <unistd.h>
#include <ncurses.h>
#include <sys/ioctl.h>
#include <stdlib.h>
#include <sys/time.h>
#include <iomanip>

//! Displays a progress bar on the screen.
/*! The ProgressBar class displays a progress bar in the console as well as an estimate of the time remaining in hour:minute:second.
 */
class ProgressBar
{
public:
    //!initializer for ProgressBar
    //! \param max maximum value corresponding to 100%
    ProgressBar(int max);
    /*! reset the progressbar such that value=0, max_value=max and start_time
        is the current time
    */
    void reset(int max);
    //!  sets current time and displays empty progress bar.
    void start();
    //!  updates progress bar with percentage value/max_value
    void update(int val);
    //!  forces progress bar to 100% and ends line
    void finish();
private:
    //!   
    void update_display();
    //!
    int value,max_value;
    //!
    double start_time,current_time;
};
#endif
