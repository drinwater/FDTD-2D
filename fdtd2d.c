#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>



int fdtd2d(const double Nx,const double Ny,double dx, double dy,double Nt,double df, double r,float l,float b, double ur,double er,double nbc,double freq,double epssrcint,double musrcint,double src){
    int NPML[4] = {20,20,20,20};
    const int Nxsize = 1000;
    const int Nysize = 1000;
    const double e0 = 8.85e-12;
    int Nx2 = 2*Nx;
    int * pNx2 = &Nx2;
    int Ny2 = 2*Ny;
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

    double ERzz;
    ERzz = (double *) malloc(sizeof(int[Nxsize][Nysize]));
    int (*URxx[Nysize]) = malloc(sizeof(int[Nxsize][Nysize]));
    int (*URyy[Nysize]) = malloc(sizeof(int[Nxsize][Nysize]));

    for (int i=0;i<Nxsize;i++){
        for (int j=0;j<Nysize;j++){
            ERzz[i][j]=0; 
        }
    }
    
    for (int i=nx1;i<nx2;i++){
        for (int j=ny1;j<ny2;j++){
            ERzz[i][j]=1; 
        }
    }

    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            if(i>=nx1 && i<=nx2 && f>=ny1 && f<=ny2){
                ERzz[i][f] = er;
            }
            else{
                ERzz[i][f] = epssrcint;
            }
        }
    }

    printf("%d\n", ERzz[nx2+1][ny1+1]);
    //UR2 = musrcint*(1-UR2) + ur*UR2;
    //Exx, Eyy, Ezz
    int o = nx1;
    FILE *fp = fopen("Edata.txt", "w");
    printf("gshgshg\n");
    if (fp == NULL)
    {
        printf("Error opening the file %s", "Edata.txt");
        return -1;
    }
    printf("File opened\n");
    // write to the text file
    for (int i = 0; i < Nxsize; i++)
    {
        for (int l = 0; l < Nysize; l++)
        {   
            
            fprintf(fp, "%d\t", ERzz[i][l]);
        }
        fprintf(fp,"\n");
    }
    printf("written\n");
    // close the file
    fclose(fp);
    printf("plotting\n");
    system("gnuplot -p -e \"set xrange[0:149];set yrange[0:149];set cbrange[-2:2];set palette defined (-2 'blue', 0 'white', 2 'red');plot 'Edata.txt' matrix with image;\"");
    getchar();
    
    //Source stuff
    printf("Source computing\n");
    double dt = 1/(Nt*freq);
    double tau = 0.5/freq;
    double t0 = 10*tau;
    int steps = 1/(dt*df);
    steps = 10000;
    double t[steps];
    for (int i = 0; i < steps; i++)
    {
        t[i] = i*dt;
    }  
    double Ezsrcint[steps];
    for (int i = 0; i < steps; i++)
    {
        Ezsrcint[i] = pow(exp(-((t[i]-t0)/tau)),2);
    }
    double nsrcint = sqrt(musrcint*epssrcint);
    int c0 = 299792458;
    double diract = (nsrcint*dy/(2*c0)) + dt/2;
    double Hxsrcint[steps];
    for (int i = 0; i < steps; i++)
    {
        Hxsrcint[i] = sqrt(epssrcint/musrcint)*pow(exp(-((t[i]+diract-t0)/tau)),2);
    }
    printf("Source computed\n");

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
    double (*sigHx[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    int e=0, g = 0;
    for (int i=0;i<Nx;i++){
        for (int f=1;f<Ny;f++){
            sigHx[i][f] = sigx[o][g];
            g+=2;
        }
        o+=2;
    }
    printf("sigHx defined");
    double (*sigHy[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    e=0, g = 0;
    for (int i=0;i<Nx;i++){
        for (int f=1;f<Ny;f++){
            sigHy[i][f] = sigy[o][g];
            g+=2;
        }
        o+=2;
    }
    double (*mHx0[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mHx0[i][f] = (1/dt) + sigHy[i][f]/(2*e0);
        }
    }
    printf("mHx0 defined");
    double (*mHx1[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mHx1[i][f] = ((1/dt) - sigHy[i][f]/(2*e0))/mHx0[i][f];
        }
    }
    printf("mHx1 defined");
    double (*mHx2[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            mHx2[i][f] = - c0/(URxx[i][f]*mHx0[i][f]);
        }
    }
    printf("mHx2 defined");
    double (*mHx3[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mHx3[i][f] = - (c0*dt/e0) * sigHx[i][f]/URxx[i][f] / mHx0[i][f];
        }
    }
    printf("mHx3 defined");
    e=0, g = 0;
    for (int i=1;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            sigHx[i][f] = sigx[o][g];
            g+=2;
        }
        o+=2;
    }
    printf("sigHx defined");
    e=0, g = 0;
    for (int i=1;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            sigHy[i][f] = sigy[o][g];
            g+=2;
        }
        o+=2;
    }
    printf("sigHy defined");
    double (*mHy0[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mHy0[i][f] = (1/dt) + sigHx[i][f]/(2*e0);
        }
    }
    printf("mHy0 defined");
    double (*mHy1[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mHy1[i][f] = ((1/dt) - sigHx[i][f]/(2*e0))/mHy0[i][f];
        }
    }
    printf("mHy1 defined");
    double (*mHy2[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    printf("Entering loop");
    for (int i=0;i<Nxsize;i++){
        printf("%d\n",i);
        for (int f=0;f<Nysize;f++){
            
            mHy2[i][f] = - c0/(URyy[i][f]*mHy0[i][f]);
        }
    }
    printf("mHy2 defined");
    double (*mHy3[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mHy3[i][f] = - (c0*dt/e0) * sigHy[i][f]/URyy[i][f] / mHy0[i][f];
        }
    }
    printf("mHy3 defined");
    double (*sigDx[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    e=0, g = 0;
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            sigDx[i][f] = sigx[o][g];
            g+=2;
        }
        o+=2;
    }
    double (*sigDy[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    e=0, g = 0;
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            sigDy[i][f] = sigy[o][g];
            g+=2;
        }
        o+=2;
    }
    double (*mDz0[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mDz0[i][f] = (1/dt) + (sigDx[i][f] + sigDy[i][f])/(2*e0) + sigDx[i][f]*sigDy[i][f]*(dt/4/(e0*e0));
        }
    }
    double (*mDz1[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mDz1[i][f] = (1/dt) - (sigDx[i][f] + sigDy[i][f])/(2*e0) - sigDx[i][f]*sigDy[i][f]*(dt/4/(e0*e0));
            mDz1[i][f] = mDz1[i][f] / mDz0[i][f];
        }
    }
    double (*mDz2[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mDz2[i][f] = c0/mDz0[i][f];
        }
    }
    double (*mDz4[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nx;i++){
        for (int f=0;f<Ny;f++){
            mDz4[i][f] = - (dt/(e0*e0))*sigDx[i][f]*sigDy[i][f]/mDz0[i][f];
        }
    }
    printf("Update coefficients computed\n");
    
    //Main FDTD loop
    printf("Defining variables");
    double (*CEx[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            CEx[i][f] = 0; 
        }
    }
    double (*CEy[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            CEy[i][f] = 0; 
        }
    }
    double (*Ez[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            Ez[i][f] = 0; 
        }
    }
    double (*Hx[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            Hx[i][f] = 0; 
        }
    }
    double (*Hy[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            Hy[i][f] = 0; 
        }
    }
    double (*CHz[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            CHz[i][f] = 0; 
        }
    }
    double (*Dz[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            Dz[i][f] = 0; 
        }
    }
    double (*ICEx[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            ICEx[i][f] = 0; 
        }
    }
    double (*ICEy[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            ICEy[i][f] = 0; 
        }
    }
    double (*IDz[Nysize]) = malloc(sizeof(double[Nxsize][Nysize]));
    for (int i=0;i<Nxsize;i++){
        for (int f=0;f<Nysize;f++){
            IDz[i][f] = 0; 
        }
    }
    printf("\nVariables defined\n");

    dx = dx*1e9;
    dy = dy*1e9;
    int T = 0;
    printf("Entering loop, T = %d\n", T);
    for(T = 0; T<=steps; T++){
        //Find Curl of Ex and Ey
        for(int nx = 1; nx<= Nx;nx++){
            for(int ny = 1;ny<=Ny-1; ny++){
                CEx[nx][ny] = (Ez[nx][ny+1] - Ez[nx][ny])/dy;
            }
            CEx[nx][Nysize] = (0 - Ez[nx][Nysize])/dy;
        }
        //Inject source to the curl of E
        for(int i = 2;i<=Nx;i++){
            CEx[i][srcint-1] = (Ez[i][srcint] - Ez[i][srcint-1])/dy - Ezsrcint[T]/dy;
        }
        for(int ny = 1; ny<=Ny;ny++){
            for(int nx = 1; nx<=Nx-1;nx++){
                CEy[nx][ny] = - (Ez[nx+1][ny] - Ez[nx][ny])/dx;
            }
        CEy[Nxsize][ny] = - (0 - Ez[Nxsize][ny])/dx;
        }
        
        //Update H integrations
        for(int nx = 1; nx<= Nx;nx++){
            for(int ny = 1;ny<=Ny; ny++){
                ICEx[nx][ny] = ICEx[nx][ny] + CEx[nx][ny];
                ICEy[nx][ny] = ICEy[nx][ny] + CEy[nx][ny];
            }
        }

        //Update H field
        for(int nx = 1; nx<= Nx;nx++){
            for(int ny = 1;ny<=Ny; ny++){
                Hx[nx][ny] = mHx1[nx][ny]*Hx[nx][ny] + mHx2[nx][ny]*CEx[nx][ny] + mHx3[nx][ny]*ICEx[nx][ny];
                Hy[nx][ny] = mHy1[nx][ny]*Hy[nx][ny] + mHy2[nx][ny]*CEy[nx][ny] + mHx3[nx][ny]*ICEy[nx][ny];
            }
        }
        //Find curl of H
        CHz[1][1] = (Hy[1][1] - 0)/dx - (Hx[1][1] - 0)/dy;
        for(int nx = 2; nx<= Nx; nx++){
            CHz[nx][1] = (Hy[nx][1] - Hy[nx-1][1])/dx - (Hx[nx][1] - 0)/dy;
        }
        for(int ny = 2; ny<= Ny; ny++){
            CHz[1][ny] = (Hy[1][ny] - 0)/dx - (Hx[1][ny] - Hx[1][ny-1])/dy;
            for(int nx = 2;nx<= Nx; nx++){
                CHz[nx][ny] = (Hy[nx][ny] - Hy[nx-1][ny])/dx - (Hx[nx][ny] - Hx[nx][ny-1])/dy;
            }
        }
        
        //Inject source to the curl of H
        for(int i = 2; i<=Nx; i++){
            CHz[i][srcint] = (Hy[i][srcint] - Hy[i-1][srcint])/dx - (Hx[i][srcint] - Hx[i][srcint-1])/dy + Hxsrcint[T]/dy;
        }
        
        //Update D integrations
        for(int nx = 1; nx<= Nx;nx++){
            for(int ny = 1;ny<=Ny; ny++){
                IDz[nx][ny] = IDz[nx][ny] + Dz[nx][ny];
            }
        }  
        
        //Update Dz
        for (int i = 0; i < Nx; i++){
            for (int l = 0; l < Ny; l++){
                Dz[i][l] = mDz1[i][l]*Dz[i][l] + mDz2[i][l]*CHz[i][l] + mDz4[i][l]*IDz[i][l];
            } 
        }
        
        //Update Ez
        for (int nx = 0; nx < Nx; nx++)
        {
            for (int ny = 0; ny < Ny; ny++)
            {
                Ez[nx][ny] = mDz1[nx][ny]*Dz[nx][ny];
            }
        }
        char filename[100] = "./data/Edata";
        char num[100];
        sprintf(num, "%d", T);
        strcat(filename,num);
        strcat(filename,".txt");

        FILE *fp = fopen(filename, "w");
        
        if (fp == NULL)
        {
            printf("Error opening the file %s", "Edata.txt");
            return -1;
        }
        
        // write to the text file
        for (int i = 0; i < Nxsize; i++)
        {
            for (int l = 0; l < Nysize; l++)
            {   
                
                fprintf(fp, "%d\t", Ez[i][l]);
            }
            fprintf(fp,"\n");
        }
        
        // close the file
        fclose(fp);
        char command[100] = "gnuplot -p -e \"load 'plotscript.gnu' ";
        strcat(command,num);
        strcat(command," \"");
        printf("T=%d\n",T);
        //system(command);
}
    return 0;

}
    
void main(){
    double NPML[] = {20, 20, 20, 20};
    double b = 0.00000005;
    float l = b;
    //int errorcode = fdtd2d(300e-9,300e-9,2e-9,2e-9,10,1e6,25e-9,l,b,-2,-2 , 1, 1e8, 1, 1,80e-9);
    printf("Starting program\n");
    getchar();
}