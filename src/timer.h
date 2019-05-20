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
