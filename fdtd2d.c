#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int fdtd2d(const double Nx,const double Ny,double dx, double dy,double Nt,double df, double r,float l,float b, double ur,double er,double nbc,double freq,double epssrcint,double musrcint,double src){
    int NPML[4] = {20,20,20,20};
    const int Nxsize = ceil(Nx/dx);
    const int Nysize = ceil(Ny/dy);
    const double e0 = 8.85e-12;
    int Nx2 = 2*Nxsize;
    int * pNx2 = &Nx2;
    int Ny2 = 2*Nysize;
    int srcint = ceil(src/dx);
    double dx2 = dx/2;
    double dy2 = dy/2;
    double epsrcint = 1;
    musrcint = 1;

    //Build geometry

    double nx = ceil(l/dx);
    int nx1 = 1 + floor((Nxsize - nx)/2);
    int nx2 = nx1 + nx - 1;
    
    double ny = ceil(b/dy);
    int ny1 = 1 + floor((Nysize - ny)/2);
    int ny2 = ny1 + ny - 1;

    //double *ERzz = (double *) malloc(Nxsize*Nysize*sizeof(double));
    size_t ERzz_length = Nxsize*Nysize;
    double *ERzz = (double *) malloc(ERzz_length*sizeof(double));
    double *URxx = (double *) malloc(Nxsize*Nysize*sizeof(double));
    double *URyy = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int j=0;j<Nysize;j++){
            if (i>=nx1 && i<=nx2 && j>=ny1 && j<=ny2)
            {
                *(ERzz + i*Nysize + j)=er;
                *(URxx + i*Nysize + j)=ur;
                *(URyy + i*Nysize + j)=ur;

            }
            else{
                *(ERzz + i*Nysize + j)=epssrcint;
                *(URxx + i*Nysize + j)=musrcint;
                *(URyy + i*Nysize + j)=musrcint;
            }
        }
    }

    //UR2 = musrcint*(1-UR2) + ur*UR2;
    //Exx, Eyy, Ezz
    FILE *fp = fopen("Edata.txt", "w");
    printf("gshgshg\n");
    if (fp == NULL)
    {
        printf("Error opening the file %s", "Edata.txt");
        return -1;
    }
    printf("File opened\n");
    // write to the text file
    for (int i= 0; i< Nxsize; i++)
    {
        for (int l = 0; l < Nysize; l++)
        {   
            
            fprintf(fp, "%f\t", ERzz[i*Nxsize + l]);
        }
        fprintf(fp,"\n");
    }
    printf("written\n");
    // close the file
    fclose(fp);
    printf("plotting\n");
    //system("gnuplot -p -e \"set xrange[0:149];set yrange[0:149];set cbrange[-2:2];set palette defined (-2 'blue', 0 'white', 2 'red');plot 'Edata.txt' matrix with image;\"");
    
    //Source stuff
    printf("Source computing\n");
    double dt = 1/(Nt*freq);
    double tau = 0.5/freq;
    double t0 = 10*tau;
    int steps = 1/(dt*df);
    steps = 10000;
    double t[steps];
    for (int index = 0; index < steps; index++)
    {
        t[index] = index*dt;
    }    
    double Ezsrcint[steps];
    for (int index = 0; index < steps; index++)
    {
        Ezsrcint[index] = exp(-pow(((t[index]-t0)/tau),2.0));
    } 
    double nsrc = sqrt(musrcint*epssrcint);
    const int c0 = 299792458;
    double diract = (nsrc*dy/(2*c0)) + dt/2;
    double Hxsrcint[steps];
    for (int index = 0; index < steps; index++)
    {
        Hxsrcint[index] =  sqrt(epssrcint/musrcint)*exp(-pow(((t[index]+diract-t0)/tau),2));
    }
    printf("Source computed\n");

    FILE *fp0 = fopen("source.txt", "w");
    printf("gshgshg\n");
    if (fp0 == NULL)
    {
        printf("Error opening the file %s", "Edata.txt");
        return -1;
    }
    printf("File opened\n");
    // write to the text file
    for (int i= 0; i< steps; i++)
    {  
            
        fprintf(fp0, "%.10f\t%.10f",  t[i],Ezsrcint[i]);
        fprintf(fp0,"\n");
    }
    printf("written\n");
    // close the file
    fclose(fp0);
    FILE *fp3 = fopen("t.txt", "w");
    printf("gshgshg\n");
    if (fp3 == NULL)
    {
        printf("Error opening the file %s", "Edata.txt");
        return -1;
    }
    printf("File opened\n");
    // write to the text file
    for (int i= 0; i< steps; i++)
    {  
            
        fprintf(fp3, "%f\t", t[i]);
    }
    printf("written\n");
    // close the file
    fclose(fp3);
    //PML
    printf("PML computing\n");
    double sigx[Nx2][Ny2];
    for (int i=0;i<Nx2;i++){
        for (int f=0;f<Ny2;f++){
            sigx[i][f] = 0; 
        }
    }
    for(int nx = 1; nx<=2*NPML[2];nx++){
        double nx1 = 2*NPML[1] - nx + 1;
        for (int i = 0; i < Ny2; i++)
        {
            sigx[nx][i] = (0.5*e0/dt)*pow((nx/2/NPML[1]),3);
        }
    }  
    for(int nx = 1; nx<=2*NPML[2];nx++){
        double nx1 = Nx2 - 2*NPML[2] + nx;
        for (int i = 0; i < Ny2; i++)
        {
            sigx[nx][i] = (0.5*e0/dt)*pow((nx/2/NPML[2]),3);
        }
    }
    double sigy[Nx2][Ny2];
    for (int i=0;i<Nx2;i++){
        for (int f=0;f<Ny2;f++){
            sigy[i][f] = 0; 
        }
    }
    for(int nx = 1; nx<=2*NPML[3];nx++){
        double ny1 = 2*NPML[3] - ny + 1;
        for (int f=0;f<Ny2;f++){
            sigy[f][nx] = (0.5*e0/dt)*pow((ny/2/NPML[3]),3);
        }
    }
    for(int nx = 1; nx<=2*NPML[4];nx++){
        double ny1 = Ny2 - 2*NPML[4] + ny;
        for (int f=0;f<Ny2;f++){
            sigy[f][nx] = (0.5*e0/dt)*pow((ny/2/NPML[4]),3);
        }
    }
    printf("PML computed\n");

    //PML update coefficients
    printf("Update coefficients computing\n");
    double *sigHx = (double *) malloc(Nxsize*Nysize*sizeof(double));
    int e=0, g=0;
    for (int i=0;i<Nxsize;i++){
        for (int f=1;f<Nysize;f++){
            *(sigHx + i*Nysize + f) = sigx[e][g];
            g+=2;
        }
        g=0;
        e+=2;
    }
    printf("sigHx defined");
    double *sigHy = (double *) malloc(Nxsize*Nysize*sizeof(double));
    e=0, g = 0;
    for (int i=0;i<Nxsize;i++){
        for (int f=1;f<Nysize;f++){
            *(sigHy + i*Nysize + f) = sigy[e][g];
            g+=2;
        }
        e+=2;
    }
    double *mHx0 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mHx0 + i*Nysize + f) = (1/dt) + *(sigHy + i*Nysize + f)/(2*e0);
        }
    }
    printf("mHx0 defined");
    double *mHx1 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mHx1 + i*Nysize + f) = ((1/dt) - *(sigHy + i*Nysize + f)/(2*e0))/ *(mHx0 + i*Nysize + f);
        }
    }
    printf("mHx1 defined\n");
    double *mHx2 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mHx2 + i*Nysize + f) = - c0/(*(URxx + i*Nysize + f)* *(mHx0 + i*Nysize + f));
        }
    }
    printf("%0.50f\n", *(mHx2 + (Nxsize-1)*Nysize + Nysize));
    printf("mHx2 defined");
    double *mHx3 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mHx3 + i*Nysize + f) = - (c0*dt/e0) * *(sigHx + i*Nysize + f)/ *(URxx + i*Nysize + f) / *(mHx0 + i*Nysize + f);
        }
    }
    printf("mHx3 defined");
    e=0, g = 0;
    for (int i=1;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(sigHx + i*Nysize + f)  = sigx[e][g];
            g+=2;
        }
        g=0;
        e+=2;
    }
    printf("sigHx defined");
    e=0, g = 0;
    for (int i=1;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(sigHy + i*Nysize + f)  = sigy[e][g];
            g+=2;
        }
        e+=2;
    }
    printf("sigHy defined");
    double *mHy0 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize; i++){
        for (int f=0;f<Nysize;f++){
            *(mHy0 + i*Nysize + f)  = (1/dt) + *(sigHx + i*Nysize + f) /(2*e0);
        }
    }
    printf("mHy0 defined");
    double *mHy1 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mHy1 + i*Nysize + f)  = ((1/dt) - *(sigHx + i*Nysize + f) /(2*e0))/ *(mHy0 + i*Nysize + f) ;
        }
    }
    
    double *mHy2 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    printf("Entering loop");
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            
            *(mHy2 + i*Nysize + f)  = - c0/(*(URyy + i*Nysize + f) * *(mHy0 + i*Nysize + f) );
        }
    }
    printf("mHy2 defined");
    double *mHy3 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mHy3 + i*Nysize + f)  = - (c0*dt/e0) * *(sigHy + i*Nysize + f) / *(URyy + i*Nysize + f)  / *(mHy0 + i*Nysize + f) ;
        }
    }
    printf("mHy3 defined");
    double *sigDx = (double *) malloc(Nxsize*Nysize*sizeof(double));
    e=0, g = 0;
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(sigDx + i*Nysize + f)  = sigx[e][g];
            g+=2;
        }
        g=0;
        e+=2;
    }
    double *sigDy = (double *) malloc(Nxsize*Nysize*sizeof(double));
    e=0, g = 0;
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(sigDy + i*Nysize + f)  = sigy[e][g];
            g+=2;
        }
        g=0;
        e+=2;
    }
    double *mDz0 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mDz0 + i*Nysize + f)  = (1/dt) + (*(sigDx + i*Nysize + f)  + *(sigDy + i*Nysize + f) )/(2*e0) + *(sigDx + i*Nysize + f) * *(sigDy+ i*Nysize + f) *(dt/4/(e0*e0));
        }
    }
    double *mDz1 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mDz1 + i*Nysize + f)  = (1/dt) - (*(sigDx + i*Nysize + f)  + *(sigDy + i*Nysize + f) )/(2*e0) - (*(sigDx + i*Nysize + f) * *(sigDy + i*Nysize + f) * (dt/4/(e0*e0)));
            *(mDz1 + i*Nysize + f)  = *(mDz1 + i*Nysize + f)  / *(mDz0 + i*Nysize + f) ;
        }
    }
    double *mDz2 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mDz2 + i*Nysize + f)  = c0/ *(mDz0 + i*Nysize + f) ;
        }
    }
    double *mDz4 = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(mDz4 + i*Nysize + f)  = - (dt/(e0*e0))* *(sigDx + i*Nysize + f) * *(sigDy + i*Nysize + f) / *(mDz0 + i*Nysize + f) ;
        }
    }
    printf("Update coefficients computed\n");
    
    //Main FDTD loop
    printf("Defining variables");
    double *CEx = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(CEx + i*Nysize + f)  = 0; 
        }
    }
    double *CEy = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(CEy + i*Nysize + f)  = 0; 
        }
    }
    double *Ez = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(Ez + i*Nysize + f)  = 0; 
        }
    }
    double *Hx = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(Hx + i*Nysize + f) = 0; 
        }
    }
    double *Hy = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(Hy + i*Nysize + f)  = 0; 
        }
    }
    double *CHz = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(CHz + i*Nysize + f)  = 0; 
        }
    }
    double *Dz = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(Dz + i*Nysize + f) = 0; 
        }
    }
    double *ICEx = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(ICEx + i*Nysize + f)  = 0; 
        }
    }
    double *ICEy = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(ICEy + i*Nysize + f)  = 0; 
        }
    }
    double *IDz = (double *) malloc(Nxsize*Nysize*sizeof(double));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            *(IDz + i*Nysize + f)  = 0; 
        }
    }
    printf("\nVariables defined\n");

    dx = dx*1e9;
    dy = dy*1e9;
    int T = 0;
    
    int index1=0, index2 = 0;
    printf("Entering loop, T = %d\n", T);
    for(T = 0; T<=steps; T++){
        //Find Curl of Ex and Ey
        for (index1 = 0; index1 < Nxsize; index1++){
            for (index2 = 0; index2< Nysize-1; index2++){
                *(CEx + index1*Nysize + index2) = (*(Ez + index1*Nysize + index2+1) - *(Ez + index1*Nysize + index2))/dy;
            }  
            *(CEx + index1*Nysize + Nysize) = (0 - *(Ez + index1*Nysize + Nysize))/dy;
        }
        //printf("%.50f\n",*(CEx + 75*Nysize + 75));
        //printf("%f\n",*(CEx + 60*Nysize + Nysize));
        //Inject source to the curl of E
        //printf("%0.50f\n",Ezsrcint[T]);
        for (index1 = 1; index1 < Nxsize; index1++){
            *(CEx + index1*Nysize + srcint-1) = (*(CEx + index1*Nysize + srcint) - *(CEx + index1*Nysize + srcint-1))/dy - Ezsrcint[T]/dy;     
        }
        
        for (index2 = 0; index2 < Nysize; index2++){
            for (index1 = 0; index1 < Nxsize-1; index1++){
                *(CEy + index1*Nysize + index2) = - (*(CEy + (index1+1)*Nysize + index2) - *(CEy + index1*Nysize + index2))/dx;
            }
            *(CEy + Nxsize + index2) = - (0 - *(CEy + Nxsize + index2))/dx;
        }
        //printf("%.50f\n",*(CEy + 75*Nysize + 75));
        //Update H integrations
        for(index1 = 0; index1<Nxsize;index1++){
            for(index2 = 0;index2<Nysize; index2++){
                *(ICEx + index1*(Nysize) + index2) = *(ICEx + index1*(Nysize) + index2) + *(CEx + index1*(Nysize) + index2);
                *(ICEy + index1*(Nysize) + index2) = *(ICEy + index1*(Nysize) + index2) + *(CEy + index1*(Nysize) + index2);
            }
        }
        //printf("%.50a\n",*(ICEx + 75*Nysize + 75));
        //Update H field
        //printf("%.50f\n",*(mHx2 + index1*(Nysize) + index2));
        for(index1 = 0; index1<Nxsize;index1++){
            for(index2 = 0;index2<Nysize; index2++){
                double xterm1 = *(mHx1 + index1*(Nysize) + index2) * (*(Hx + index1*(Nysize) + index2));
                double xterm2 = *(mHx2 + index1*(Nysize) + index2) * (*(CEx + index1*(Nysize) + index2));
                double xterm3 = *(mHx3 + index1*(Nysize) + index2) * (*(ICEx + index1*(Nysize) + index2));
                *(Hx + index1*(Nysize) + index2) =  xterm1 + xterm2 + xterm3 ;
                double hxterm = *(Hx + index1*(Nysize) + index2);
                *(Hy + index1*(Nysize) + index2) = *(mHy1 + index1*(Nysize) + index2) * *(Hy + index1*(Nysize) + index2) + *(mHy2 + index1*(Nysize) + index2) * *(CEy + index1*(Nysize) + index2) + *(mHy3 + index1*(Nysize) + index2) * *(ICEy + index1*(Nysize) + index2);
            }
        }
        //printf("%.50a\n",*(Hx + 75*Nysize + 75));
        //Find curl of H
        *(CHz + 1*Nysize + 1) = (*(Hy + 1*Nysize + 1) - 0)/dx - (*(Hx + 1*Nysize + 1) - 0)/dy;
        for(index1 = 1; index1< Nxsize; index1++){
            *(CHz + index1*Nysize + 1) = (*(Hy + index1*Nysize + 1) - *(Hy + (index1-1)*Nysize + 1))/dx - (*(Hx + index1*Nysize + 1) - 0)/dy;
        }
            
        for(index2 = 1; index2<Nysize; index2++){
            *(CEx + Nysize + index2) = (*(Hy + Nysize + index2) - 0)/dx - (*(Hx + Nysize + index2) - *(Hx + Nysize + index2-1))/dy;
            for(index1 = 1; index1<Nxsize; index1++){
                *(CHz + index1*Nysize + index2) = (*(Hy + index1*Nysize + index2) - *(Hy + (index1-1)*Nysize + index2))/dx - (*(Hx + index1*Nysize + index2) - *(Hx + index1*Nysize + index2-1))/dy;
            }  
        }
        //printf("%.50f\n",*(CHz + 75*Nysize + 75));
        //Inject source to the curl of H
        for(index1 = 1; index1<Nxsize; index1++){
            *(CHz + index1*Nysize + srcint) = (*(Hy + index1*Nysize + srcint) - *(Hy + (index1-1)*Nysize + srcint))/dx - (*(Hx + index1*Nysize + srcint) - *(Hx + index1*Nysize + srcint-1))/dy + Hxsrcint[T]/dy;
        }
        //Update D integrations
        for(index1 = 0; index1< Nxsize;index1++){
            for(index2 = 0;index2<Nysize; index2++){
                *(IDz + index1*Nysize + index2) = *(IDz + index1*Nysize + index2) + *(Dz + index1*Nysize + index2);
            }
        }  
        //printf("%.50f\n",*(IDz + 75*Nysize + 75));
        //Update Dz
        for (index1 = 0; index1 < Nxsize; index1++){
            for (index2 = 0; index2 < Nysize; index2++){
                *(Dz + index1*Nysize + index2) = *(mDz1 + index1*Nysize + index2) * *(Dz + index1*Nysize + index2) + *(mDz2 + index1*Nysize + index2) * *(CHz + index1*Nysize + index2) + *(mDz4 + index1*Nysize + index2) * *(IDz + index1*Nysize + index2);
            } 
        }
        //printf("%.50f",*(Dz + 75*Nysize + 75));
        //Update Ez
        for (index1 = 0; index1 < Nxsize; index1++)
        {
            for (index2 = 0; index2 < Nysize; index2++)
            {
                *(Ez + index1*Nysize + index2) = *(mDz1 + index1*Nysize + index2) * *(Dz + index1*Nysize + index2);
            }
        }
        //printf("%.50f\n",*(Ez + 75*Nysize + 75));
        char filename[100] = "./data/Edata";
        char num[100];
        sprintf(num, "%d", T);
        strcat(filename,num);
        strcat(filename,".txt");

        FILE *fp2 = fopen(filename, "w");
        
        if (fp2 == NULL)
        {
            printf("Error opening the file %s", "Edata.txt");
            return -1;
        }
        
        // write to the text file
        for (index1= 0; index1< Nxsize; index1++)
        {
            for (index2 = 0; index2 < Nysize; index2++)
            {   
                fprintf(fp2, "%.50f\t", *(Ez + index1*Nysize + index2));
            }
            fprintf(fp2,"\n");
        }
        
        // close the file
        fclose(fp2);
        char command[100] = "gnuplot -p -e \"load 'plotscript.gnu' ";
        strcat(command,num);
        strcat(command," \"");
        //system(command);
    }
    return 0;
}
    
void main(){
    double NPML[] = {20, 20, 20, 20};
    double b = 0.00000005;
    float l = b;
    printf("Starting program\n");
    //int errorcode = fdtd2d(300e-9,300e-9,2e-9,2e-9,10,1e6,25e-9,l,b,-2,-2 , 1, 1e8, 1, 1,80e-9);
    //printf("Errorcode:%d",errorcode);
    //getchar();
}