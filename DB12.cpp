#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <iostream>
#include <string.h>
#include <numeric>
#include <algorithm>
#include <string>

const double NV_MAGICCONST = 4 * exp(-0.5)/sqrt(2.0);

bool DEBUG = false;

double sampleNormal() {
  double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
  double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
  double r = u * u + v * v;
  if (r == 0 || r > 1) return sampleNormal();
  double c = sqrt(-2 * log(r) / r);
  return u * c;
}

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

float singleDiracBenchmark(pid_t pid){

  int initialIterations=2;
  int iterations = 1;
  int extraIterations = 2;
  //  std::cout << "Function running in thread " << thread_id << std::endl;

  int nTOT = 1000 * 1000 * 12.5 ;
  float calib = 25.0;
  clock_t t;
  float cput;

  timespec time1, time2;

  for (int i=0; i< initialIterations + iterations + extraIterations; ++i){
    //std::cout << std::endl << "i " << i << std::endl;

    int first_time;
    if (i==initialIterations){
      if(DEBUG) std::cout << pid << " Starting the measurement loop..." << std::endl;
      first_time = time(NULL);
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time1);
    }
    
    for(int n=0; n<nTOT; ++n) {
      //sampleNormal();
      normalvariate(); //This implements the same function used in the python version
    }

    if (i == initialIterations+iterations-1){
      int second_time = time(NULL);
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time2);
      cput = diff(time1,time2);

      if(DEBUG) std::cout << pid << " Ending the measurement loop... It took me " << second_time - first_time << " second and cput " << cput << std::endl;
    }
  }

  if (cput==0)
    return 99999999999;
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


  pid_t childpid;
  int my_process_id = 0;
  int pipe_file_desc[2];
  int fd[2*num_of_threads];

  if(pipe(pipe_file_desc) == -1)
    {perror("error creating pipe");
      return 1;
    }

  

  for (int i=0;i<num_of_threads;++i){

    if(pipe(&fd[2*i]) == -1)
      {perror("error creating pipe");
	return 1;
      }

    if((childpid = fork()) == -1)
      {
	perror("fork() failed\n");
	exit(1);
      }    

    if (childpid == 0){
      //Child

      read(pipe_file_desc[0], &my_process_id,4);
      if(DEBUG) printf("child process %d ID: %d\n",i,my_process_id);

      // child process
      float db12value = singleDiracBenchmark(my_process_id);

      char astring[120] ;
      sprintf(astring,"%f",db12value);
      //std::cout << my_process_id << " db12 " << astring << std::endl;
      write(fd[2*i+1], astring, (strlen(astring)+1));

      exit(0);
    }else{ //childpid>0
      // parent process
      if(DEBUG) printf("parent has child %d with pid %d\n",i,childpid);
      write(pipe_file_desc[1], &childpid,4);
    }
  }

  char readbuffer[80];
  std::vector<float> results;
  for(int i=0; i<num_of_threads; i++){
    read(fd[2*i], readbuffer, sizeof(readbuffer));
    results.push_back(atof(readbuffer));
    if(DEBUG) printf("Parent Received string: %s\n", readbuffer);
  }

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
