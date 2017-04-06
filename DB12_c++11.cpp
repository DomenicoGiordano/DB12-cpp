#include <unistd.h>
#include <iostream>
#include <random>
#include <cmath>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <thread>
#include <future>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cstring>


std::random_device rd;
std::mt19937 gen(1234567890);

const double NV_MAGICCONST = 4 * exp(-0.5)/sqrt(2.0);

bool DEBUG = false;

float diff(timespec start, timespec end)
{
  timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp.tv_sec+temp.tv_nsec/1000000000.;
}

float normalvariate(double  mu=0, double sigma=1){
  /*Normal distribution.

        mu is the mean, and sigma is the standard deviation.

  */
  double u1,u2,z,zz;
  while(1){
    u1 = ((double) rand() / (RAND_MAX)) ;
    u2 = 1 - ((double) rand() / (RAND_MAX)) ;
    
    z = NV_MAGICCONST*(u1-0.5)/u2;
    zz = z*z/4.0;

    if (zz <= -log(u2))
      break;
  }
  return mu + z*sigma;
}

float singleDiracBenchmark(std::promise<float> * promObj){

  std::thread::id thread_id = std::this_thread::get_id();

  int initialIterations=2;
  int iterations = 1;
  int extraIterations = 1;
  if (DEBUG) std::cout << "Function running in thread " << thread_id << std::endl;

  int nTOT = 1000 * 1000 * 12.5 ;
  float calib = 250.0;
  clock_t t;
  float cput;

  timespec time1, time2;

  for (int i=0; i< initialIterations + iterations + extraIterations; ++i){
    //std::cout << std::endl << "i " << i << std::endl;

    int first_time;
    if (i==initialIterations){
      if (DEBUG) std::cout << "thread " << thread_id << " Starting the measurement loop..." << std::endl;
      t = clock();
      first_time = time(NULL);

       clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time1);
    }
    
   
    std::normal_distribution<> d(10,1);
    for(int n=0; n<nTOT; ++n) {
      d(gen);
    }
    
    /*
    for(int n=0; n<nTOT; ++n) {
      //sampleNormal();
      normalvariate(); //This implements the same function used in the python version
    }
    */
    if (i == initialIterations+iterations-1){
      t = clock() - t;
      int second_time = time(NULL);
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time2);
      cput = diff(time1,time2);

      if (DEBUG)  std::cout << "Ending the measurement loop... It took me " << t << " clicks (" << ((float)t)/CLOCKS_PER_SEC << " seconds on " << CLOCKS_PER_SEC << " CLOCKS_PER_SEC). and actual time " << second_time - first_time << " and cput " << cput << std::endl;
    }
  }

  promObj->set_value(calib*iterations/cput);
  return calib*iterations/cput;
}

int main(int argc, char **argv)
{
  printf("--beginning of DB12 cpp version\n");
  
  int c;
  char *Nvalue=NULL;

  opterr = 0;

  while ((c = getopt (argc, argv, "dn:")) != -1)
    switch (c)
      {
      case 'd':
        DEBUG = true;
        break;
      case 'n':
        Nvalue = optarg;
        break;
      case '?':
        if (optopt == 'n')
          fprintf (stderr, "Option -%n requires as argument the number of threads, or wholenode \n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Syntax: %s [-n <num_threads> (wholenode)] [-d (DEBUG)]\n", argv[0]);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

  int num_of_threads = sysconf(_SC_NPROCESSORS_ONLN);
  if(Nvalue!=NULL){
    short nthreads = atoi(Nvalue);
    if (nthreads > 0 and nthreads <num_of_threads)
      num_of_threads = nthreads;
  }
  std::cout << num_of_threads << " concurrent threads are supported." << std::endl;

  std::thread tt[num_of_threads];
  std::promise<float> promiseObj[num_of_threads];
  std::future<float> futureObj[num_of_threads];

  //This statement will launch multiple threads in loop
  for (int i = 0; i < num_of_threads; ++i) {
    futureObj[i]= promiseObj[i].get_future();
    tt[i] = std::thread(singleDiracBenchmark,&(promiseObj[i]));
  }

  //This is how we join the new thread with main
  for (int i = 0; i < num_of_threads; ++i) {
    tt[i].join();
  }

  std::vector<float> results;
  for (int i = 0; i < num_of_threads; ++i) {
    results.push_back(futureObj[i].get());
  }
  std::cout <<std::endl;

  float average = std::accumulate(results.begin(),results.end(), 0.0)/ results.size();
  std::sort(results.begin(), results.end());
  float min = results[0];
  float max = results.back();
  float median = results.size()%2 ==0 ? (results[results.size()/2]+results[results.size()/2-1])/2. : results[results.size()/2];
  
  char cres[1024];
  sprintf(cres, "%.2f,",results[0]);
  for (int i=1; i<results.size(); ++i)
    sprintf(strchr(cres, '\0'), "%.2f,",results[i]);
  printf("JSON: {\"ncores\": %d, \"max\": %.2f, \"min\": %.2f, \"mean\": %.2f, \"median\": %.2f, \"values\": [%s]}",num_of_threads,max,min,average,median,cres);

  
  printf("\n--end of program--\n");

  return 0;
}
