/*
Fractal code for CS 4380 / CS 5351

Copyright (c) 2018, Texas State University. All rights reserved.

Redistribution and usage in source and binary form, with or without
modification, is only permitted for educational use.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Martin Burtscher
*/

#include <cstdlib>
#include <sys/time.h>
#include "fractalpar.h"
#include <iostream>
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#define main_process 0
static const double Delta = 0.001;
static const double xMid =  0.23701;
static const double yMid =  0.521;

int main(int argc, char *argv[])
{
  MPI_Status status;
  int num_processes,process_rank, source, destination, tag;
  
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
 
  // check command line
  if (argc != 3) {if(process_rank==0)fprintf(stderr, "usage: %s frame_width num_frames\n", argv[0]); exit(-1);}
  int width = atoi(argv[1]);
  if (width < 10) {if(process_rank==0)fprintf(stderr, "error: frame_width must be at least 10\n"); exit(-1);}
  int frames = atoi(argv[2]);
  if (frames < 1) {if(process_rank==0)fprintf(stderr, "error: num_frames must be at least 1\n"); exit(-1);}
  
  if(process_rank==0) printf("Fractal v1.6 [serial|parallel]\n");
  if(process_rank==0) printf("computing %d frames of %d by %d fractal\n", frames, width, width);

  // allocate picture array
  unsigned char* pic = new unsigned char[frames * width * width];

  // start time
  timeval start, end;
  gettimeofday(&start, NULL);


  // compute frames
  double delta = Delta;
  int fatia = frames/(num_processes-1);
  int restante = frames - ((num_processes-1) * fatia);
  if(process_rank != 0) {
    restante = num_processes - process_rank != 1? 0 : restante;
    std::cout << "rank: " << process_rank << " restante : " << restante << std::endl; 
    for (int frame = 0 + (process_rank-1)*(fatia); frame < (process_rank*fatia)+restante; frame++) {
      delta = Delta * pow(0.98,frame);
      std::cout << "Frame : " << frame << "| rank : " << process_rank << "| fatia : " << fatia+restante << std::endl;
      const double xMin = xMid - delta;
      const double yMin = yMid - delta;
      const double dw = 2.0 * delta / width;
      for (int row = 0; row < width; row++) {
        const double cy = yMin + row * dw;
        for (int col = 0; col < width; col++) {
          const double cx = xMin + col * dw;
          double x = cx;
          double y = cy;
          int depth = 256;
          double x2, y2;
          do {
            x2 = x * x;
            y2 = y * y;
            y = 2 * x * y + cy;
            x = x2 - y2 + cx;
            depth--;
          } while ((depth > 0) && ((x2 + y2) < 5.0));
          pic[frame * width * width + row * width + col] = (unsigned char)depth;
        }
      }
    }
    MPI_Send(&pic[(process_rank-1)*(fatia) * width * width],(fatia+restante)*width*width,MPI_UNSIGNED_CHAR,0,0,MPI_COMM_WORLD);
  }
  else
  {
    int restanteaux = restante;
    for (int source = 1; source < num_processes; source++)
    {
      restante = num_processes-source != 1 ? 0 : restanteaux; 
      MPI_Recv(&pic[(source-1)*(fatia) * width * width],(fatia+restante)*width*width,MPI_UNSIGNED_CHAR,source,0,MPI_COMM_WORLD,&status);
    }
  }
  
  // end time
  if(process_rank==0){
    gettimeofday(&end, NULL);
    double runtime = end.tv_sec + end.tv_usec / 1000000.0 - start.tv_sec - start.tv_usec / 1000000.0;
    printf("compute time: %.4f s\n", runtime);
    if ((width <= 256) && (frames <= 100)) {
      std::cout<<"Gravando imagens"<<std::endl;
      for (int frame = 0; frame < frames; frame++) {
        char name[32];
        sprintf(name, "fractal%d.bmp", frame + 1000);
        writeBMP(width, width, &pic[frame * width * width], name);
      }
    }
  }
  // verify result by writing frames to BMP files

  delete [] pic;
  MPI_Finalize();
  return 0;
}
