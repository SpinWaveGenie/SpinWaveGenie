/* main.cpp 
* Author : Hayi NN 
* YK, 16/03/11 
*/ 
#include "progressbar.h"

using namespace std;

ProgressBar::ProgressBar(int n)
{
    max_value = n;
}

void ProgressBar::start()
{
    struct timeval tim;
    gettimeofday(&tim, NULL);
    start_time=tim.tv_sec+(tim.tv_usec/1000000.0);
    value = 0;
    update_display();
}

void ProgressBar::update(int i)
{
    struct timeval tim;
    gettimeofday(&tim, NULL);
    current_time=tim.tv_sec+(tim.tv_usec/1000000.0);
    value = i;
    update_display();
}

void ProgressBar::finish()
{
    value = max_value;
    update_display();
    cout << endl;
}

void ProgressBar::update_display()
{
    double percentage = 100.0*value/max_value;
    int pct = (int) percentage;
    string bar; 
    struct winsize uk; 
    if(ioctl(0,TIOCGWINSZ,&uk) != 0) { 
        exit(1); 
    } 
    int wdt = uk.ws_col - 30; // create some space for the other text
    if(wdt < 5) { 
        wdt = 5;      // minimum width for progress bar 5 char 
    } 
        
    for(int i=0; i < wdt; i++) { 
        if(i < (pct*wdt/100)) { 
            bar.replace(i,1,"="); 
        } else if(i == (pct*wdt/100)) { 
            bar.replace(i,1,">"); 
        } else { 
            bar.replace(i,1," "); 
        } 
    } 
    cout << "\r";          // go to the first character in terminal 
    cout << setw(3) << pct << "%";
    cout << "[" << bar << "] ";     // progress bar
       
    // add ETA to progress bar
    cout << "ETA: ";

    // get current time and find ellapsed hours, minutes, seconds
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double current_time = tim.tv_sec+(tim.tv_usec/1000000.0);
    double secperstep = (current_time-start_time)/percentage;
    double remaining_time = secperstep*(100.0-percentage);
    int ti = (int) remaining_time;
    int hour = ti/3600; ti = ti%3600;
    if (hour < 0)
        hour = 0;
    int min = ti/60; ti = ti%60;
    if (min < 0)
        min = 0;
    int sec = ti;
    if (sec < 0)
        sec = 0;
    //print out time
    cout << setw(2) << setfill('0') << hour << ":"
         << setw(2) <<  min << ":"
         << setw(2) <<  sec << flush;
    cout << setfill(' ');
}
     
/*int main(int argv, char* argc[]) {
    int n=10;
    ProgressBar pbar(n);
    pbar.start();
    for(int i =0; i < n; i++) {
        pbar.update(i);
        sleep(1);
    }
    pbar.finish();
    return 0; 
}*/ 
