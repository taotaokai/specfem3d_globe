#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))

typedef float realw;

__global__ void compute_teleseismic_gradient_kernel(realw* accel,
                                                    int interface_type,
                                                    int num_teleseismic_boundary_faces,
                                                    int* teleseismic_boundary_ispec,
                                                    realw* teleseismic_boundary_area,
                                                    int* ibool,
                                                    realw* d_field_teleseismic,
                                                    realw scale_displ) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;
  realw area;

  // don't compute surface faces outside of range
  // and don't compute points outside NGLLSQUARE==NGLL2==25
  //if(igll < NGLL2 && iface < num_teleseismic_boundary_faces) {

  // way 2: only check face, no further check needed since blocksize = 25
  if (iface < num_teleseismic_boundary_faces){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = teleseismic_boundary_ispec[iface]-1;

    // determines indices i,j,k depending on teleseismic boundary type
    switch( interface_type ){
      case 0:
        // xmin
        i = 0; // index -1
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);
        break;
      case 1:
        // xmax
        i = NGLLX-1;
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);
        break;
      case 2:
        // ymin
        j = 0;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);
        break;
      case 3:
        // ymax
        j = NGLLX-1;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);
        break;
      case 4:
        // zmin
        k = 0;
        j = (igll/NGLLX);
        i = (igll-j*NGLLX);
        break;
    }

    iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    // gets associated area
    area = teleseismic_boundary_area[INDEX2(NGLL2,igll,iface)];

    d_field_teleseismic[INDEX3(NDIM,NGLL2,0,igll,iface)] = &
      accel[iglob*3]*scale_displ*area;

    d_field_teleseismic[INDEX3(NDIM,NGLL2,1,igll,iface)] = &
      accel[iglob*3+1]*scale_displ*area;

    d_field_teleseismic[INDEX3(NDIM,NGLL2,2,igll,iface)] = &
      accel[iglob*3+2]*scale_displ*area;

  } // num_teleseismic_boundary_faces
}

