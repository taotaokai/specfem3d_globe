module BOAST
  def BOAST::compute_seismograms_kernel( ref = true, n_dim = 3, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_seismograms_kernel"

    nrec_local             = Int("nrec_local",             :dir => :in)
    displ                  = Real("displ",                 :dir => :in,   :dim => [ Dim() ])
    d_ibool                = Int("d_ibool",                :dir => :in,   :dim => [ Dim() ])
    hxir                   = Real("hxir",                  :dir => :in,   :dim => [ Dim() ])
    hetar                  = Real("hetar",                 :dir => :in,   :dim => [ Dim() ])
    hgammar                = Real("hgammar",               :dir => :in,   :dim => [ Dim() ])
    seismograms            = Real("seismograms",           :dir => :inout,:dim => [ Dim() ])
    nu                     = Real("nu",                    :dir => :in,   :dim => [ Dim() ])
    ispec_selected_rec     = Int("ispec_selected_rec",     :dir => :in,   :dim => [ Dim() ])
    number_receiver_global = Int("number_receiver_global", :dir => :in,   :dim => [ Dim() ])
    scale_displ            = Real("scale_displ",           :dir => :in)

    ndim         = Int("NDIM",         :const => n_dim)
    ngllx        = Int("NGLLX",        :const => n_gllx)
    ngll2        = Int("NGLL2",        :const => n_gll2)
    ngll3        = Int("NGLL3",        :const => n_gll3)
    ngll3_padded = Int("NGLL3_PADDED", :const => n_gll3_padded)

    p = Procedure(function_name, [nrec_local,displ,d_ibool,hxir,hetar,hgammar,seismograms,nu,ispec_selected_rec,number_receiver_global,scale_displ])

    if (get_lang == CUDA and ref) then
      get_output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx, :ngll2 => n_gll2, :ngll3 => n_gll3, :ngll3_padded => n_gll3_padded )
      open p
      ispec =      Int("ispec")
      iglob =      Int("iglob")
      irec_local = Int("irec_local")
      irec =       Int("irec")
      tx =         Int("tx")
      lagrange=    Real("lagrange")
      i =          Int("i")
      j =          Int("j")
      k =          Int("k")
      l =          Int("l")
      s =          Int("s")
      decl ispec
      decl iglob
      decl irec_local
      decl irec
      decl tx
      decl lagrange
      decl i
      decl j
      decl k
      decl l
      decl s

      decl sh_dxd = Real("sh_dxd",     :local => true, :dim => [Dim(ngll3_padded)] )
      decl sh_dyd = Real("sh_dyd",     :local => true, :dim => [Dim(ngll3_padded)] )
      decl sh_dzd = Real("sh_dzd",     :local => true, :dim => [Dim(ngll3_padded)] )

      print tx === get_local_id(0)
      print irec_local === get_group_id(0) + get_num_groups(0)*get_group_id(1)

      print k === tx/ngll2
      print j === (tx - k*ngll2)/ngllx
      print i ===  tx - k*ngll2 - j*ngllx

      print If (irec_local < nrec_local) {

        print irec  === number_receiver_global[irec_local] - 1
        print ispec === ispec_selected_rec[irec] - 1

        print sh_dxd[tx] === 0
        print sh_dyd[tx] === 0
        print sh_dzd[tx] === 0

        print If (tx < ngll3) {
          print lagrange === hxir[INDEX2(ngllx,i,irec_local)] * hetar[INDEX2(ngllx,j,irec_local)] * hgammar[INDEX2(ngllx,k,irec_local)]
          print iglob === d_ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)]-1

          print sh_dxd[tx] === lagrange * displ[INDEX2(ndim,0,iglob)]
          print sh_dyd[tx] === lagrange * displ[INDEX2(ndim,1,iglob)]
          print sh_dzd[tx] === lagrange * displ[INDEX2(ndim,2,iglob)]
        }
        print barrier(:local)

        # sum(sh_dxd[:]) is calculated in cascade of adding 2 numbers at each step
        # 0,1->0, 2,3->2, 4,5->4, 6,7->6, ...
        # 0,      2->0,   4,      6->4, ... 
        # 0,              4->0, ...
        # ...
        # here after 7 steps, indices from 0 to 127 (2**7=128 numbers) are summed into index 0
        # this is enough for ngll[x|y|z]=5, that is 5**3=125 numbers. 
        print l === 1
        (1..7).each { |indx1|
          print s === l*2
          print If ( BOAST::Modulo(tx,s) == 0) {
            print sh_dxd[tx] === sh_dxd[tx] + sh_dxd[tx + l]
            print sh_dyd[tx] === sh_dyd[tx] + sh_dyd[tx + l]
            print sh_dzd[tx] === sh_dzd[tx] + sh_dzd[tx + l]
          }
          print barrier(:local)
          print l ===  l*2
        }

        (0..2).each { |indx|
          print If (tx == indx) {
            print seismograms[INDEX2(ndim,indx,irec_local)] === scale_displ * ( nu[INDEX3(ndim,ndim,indx,0,irec_local)]*sh_dxd[0] + nu[INDEX3(ndim,ndim,indx,1,irec_local)]*sh_dyd[0] + nu[INDEX3(ndim,ndim,indx,2,irec_local)]*sh_dzd[0] )
          }
        }
        #print If (tx == 0) {
        #  print seismograms[INDEX2(ndim,0,irec_local)] === scale_displ * ( nu[INDEX3(ndim,ndim,0,0,irec_local)]*sh_dxd[0] + nu[INDEX3(ndim,ndim,0,1,irec_local)]*sh_dyd[0] + nu[INDEX3(ndim,ndim,0,2,irec_local)]*sh_dzd[0] )
        #}
        #print If (tx == 1) {
        #  print seismograms[INDEX2(ndim,1,irec_local)] === scale_displ * ( nu[INDEX3(ndim,ndim,1,0,irec_local)]*sh_dxd[0] + nu[INDEX3(ndim,ndim,1,1,irec_local)]*sh_dyd[0] + nu[INDEX3(ndim,ndim,1,2,irec_local)]*sh_dzd[0] )
        #}
        #print If (tx == 2) {
        #  print seismograms[INDEX2(ndim,2,irec_local)] === scale_displ * ( nu[INDEX3(ndim,ndim,2,0,irec_local)]*sh_dxd[0] + nu[INDEX3(ndim,ndim,2,1,irec_local)]*sh_dyd[0] + nu[INDEX3(ndim,ndim,2,2,irec_local)]*sh_dzd[0] )
        #}
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env(:array_start)
    kernel.procedure = p
    return kernel
  end

end
