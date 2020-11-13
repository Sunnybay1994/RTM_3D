#include "fdtd.h"
using namespace std;

void output_gather()
{
    int it;
    FILE *fp;
    char filename[30];
    char dir[MAX_NAME_LEN] = "./Output/";
    sprintf(filename, "gather_%02d.dat", myRank);
    fp=fopen(strcat(strcat(dir, filename),tag), "w+");
    if(fp==NULL)
    {
        printf("File open error in output_gather!");
	exit(1);
    }
    for(ii=0; ii<(int) myrec.size(); ++ii)
    {
        i = rec[ myrec[ii] ].x;
        j = rec[ myrec[ii] ].y;
        k = rec[ myrec[ii] ].z;
        fprintf(fp,"%d %d %d ",i,j,k);
        for(it=0; it<nt; ++it)
        {
            fprintf(fp, "%e ", gather(ii,it));

        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
void output_slice(int it,double *EH_total)
{
    int data_size = sizeof(EH_slx[0]);
    

    char xfilename[30],yfilename[30],zfilename[30];
    char xdir[MAX_NAME_LEN] = "./Output/";
    char ydir[MAX_NAME_LEN] = "./Output/";
    char zdir[MAX_NAME_LEN] = "./Output/";
    
    sprintf(xfilename, "slx_Ey%s_%05d.bin", tag,it);
    sprintf(yfilename, "sly_Ey%s_%05d.bin", tag,it);
    sprintf(zfilename, "slz_Ey%s_%05d.bin", tag,it);
    FILE *fp;
    
    fp=fopen(strcat(xdir, xfilename), "wb+");
    if(fp==NULL)
    {
        printf("File open error in output_slice!");
	exit(1);
    }
    if(nxslice != 0){
    for(i=0;i<nxslice;++i){
        for(k=0;k<nz;++k){
            for(j=0;j<ny;++j){
                EH_slx[(i*nz+k)*ny+j] = EH_total(slicex[i]+order,j,k);
            //     fprintf(fp, "%e ", EH_total(slicex[i]+order,j,k));
            }
            }}}
    fwrite(EH_slx,data_size,ny*nz,fp);
    fclose(fp);
    fp=fopen(strcat(ydir, yfilename), "wb+");
    if(fp==NULL)
    {
        printf("File open error");
    }
    if(nyslice !=0){
    for(j=0;j<nyslice;++j){
       for(k=0;k<nz;++k){
            for(i=0;i<nx;++i){
                EH_sly[(j*nz+k)*nx+i] = EH_total(i+order,slicey[j],k);
                // fprintf(fp, "%e ", EH_total(i+order,slicey[j],k));
            }
            }}}
    fwrite(EH_sly,data_size,nx*nz,fp);
    fclose(fp);

    fp=fopen(strcat(zdir, zfilename), "wb+");
    if(fp==NULL)
    {
        printf("File open error");
    }
    if(nzslice != 0){
    for(k=0;k<nzslice;++k){
        for(j=0;j<ny;++j){
            for(i=0;i<nx;++i){
                EH_slz[(k*ny+j)*nx+i] = EH_total(i+order,j,slicez[k]);
                // fprintf(fp, "%e ", EH_total(i+order,j,slicez[k]));
            }}}}
    fwrite(EH_slz,data_size,nx*ny,fp);
    fclose(fp);

}
void output_wavefield(int it,double *EH_total)
{
    int data_size = sizeof(EH_wvf[0]);

    char filename[30];
    char dir[MAX_NAME_LEN] = "./Output/";
    sprintf(filename, "wvf_Ey%s_%05d.bin", tag, it);
    FILE *fp;
    fp=fopen(strcat(dir,filename), "wb+");

    if(fp==NULL)
    {
        printf("File open error in output_wavefield!");
	    exit(1);
    }
    // store element in the first dimention first, to coincide with fortran
    int inz = int(nz/output_step_x_of_wavefield);
    int iny = int(ny/output_step_x_of_wavefield);
    int inx = int(nx/output_step_x_of_wavefield);
    for(k=0;k<nz;k+=output_step_x_of_wavefield)
    {
        int ik = int(k/output_step_x_of_wavefield);
        for(j=0;j<ny;j+=output_step_x_of_wavefield)
        {
            int ij = int(j/output_step_x_of_wavefield);
            for(i=0;i<nx;i+=output_step_x_of_wavefield)
            {
                int ii = int(i/output_step_x_of_wavefield);
                EH_wvf[((ik*inz+ij)*iny+ii)] = EH_total(i+order,j,k);
                // fprintf(fp, "%e ", EH_total(i+order,j,k));
            }
        }
    }
    fwrite(EH_wvf,data_size,inx*iny*inz,fp);
    fclose(fp);
}

