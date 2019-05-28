//************************************************//
// Copyright 2016-2019 University of Warwick

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
// sell copies of the Software, and to permit persons to whom the Software is furnished 
// to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//************************************************//


#ifndef TIMER_H
#define TIMER_H

typedef std::map<std::string, std::vector<double> > t_times;
t_times times;

void init_timers(std::string timers[], int length)
{
    for(int i = 0; i < length; i++)
    {
        times[timers[i]].resize(3);
    }
}

void start_timer(const std::string& timer_name, double time)
{
    times[timer_name][0] = time;
}

void end_timer(const std::string& timer_name, double time)
{
    std::vector<double>& tuple = times[timer_name];
    tuple[1] += time - tuple[0];
    tuple[2]++;
}

void dump_timers()
{
    t_times::iterator iterator;

    op_printf("Kernel runtimes:\n");
    for(iterator = times.begin(); iterator != times.end(); ++iterator)
    {
        // std::cout << iterator->first  << ": " << (iterator->second[1])
        //     << std::endl;
        op_printf("%s: %f\n", iterator->first.c_str(), iterator->second[1]);
    }
}

#endif 
