module BOAST
  def BOAST::compute_add_teleseismic_sources_kernel(ref = false, n_dim = 3, n_gllx = 5, n_gll2 = 25)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "compute_add_teleseismic_sources_kernel"

    accel = Real("accel", :dir => :inout, :dim => [ Dim() ])

    field_teleseismic = Real("field_teleseismic", :dir => :in, :dim => [ Dim() ])

    interface_type = Int("interface_type", :dir => :in)

    num_teleseismic_boundary_faces = Int("num_teleseismic_boundary_faces", :dir => :in)

    teleseismic_boundary_ispec = Int("teleseismic_boundary_ispec", :dir => :in, :dim => [ Dim() ])

    ibool = Int( "ibool", :dir => :in, :dim => [ Dim() ])

    variables = [accel, field_teleseismic, interface_type, num_teleseismic_boundary_faces, teleseismic_boundary_ispec, ibool]

    ndim  = Int("NDIM",  :const => n_dim)
    ngllx = Int("NGLLX", :const => n_gllx)
    ngll2 = Int("NGLL2", :const => n_gll2)

    p = Procedure(function_name, variables)

    if (get_lang == CUDA and ref) then

      get_output.print File::read("references/#{function_name}.cu")

    elsif(get_lang == CL or get_lang == CUDA) then

      make_specfem3d_header( :ndim => n_dim, :ngllx => n_gllx, :ngll2 => n_gll2 )
      open p
      decl igll = Int("igll")
      decl iface = Int("iface")
      decl i = Int("i"), j = Int("j"), k = Int("k")
      decl iglob = Int("iglob")
      decl ispec = Int("ispec")

      print igll === get_local_id(0)
      print iface === get_group_id(0)+get_group_id(1)*get_num_groups(0)

      print If(iface < num_teleseismic_boundary_faces ) {
        print ispec === teleseismic_boundary_ispec[iface]-1

        print Case( interface_type,
        0 => lambda {
          print i === 0
          print k === igll/ngllx
          print j === igll-k*ngllx
        },
        1 => lambda {
          print i === ngllx-1
          print k === igll/ngllx
          print j === igll-k*ngllx
        },
        2 => lambda {
          print j === 0
          print k === igll/ngllx
          print i === igll-k*ngllx
        },
        3 => lambda {
          print j === ngllx-1
          print k === igll/ngllx
          print i === igll-k*ngllx
        },
        4 => lambda {
          print k === 0
          print j === igll/ngllx
          print i === igll-j*ngllx
        })

        print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1

        (0..2).each { |indx|
          print atomicAdd(accel + iglob*3 + indx, field_teleseismic[INDEX3(ndim,ngll2,indx,igll,iface)])
        }

      }

      close p

    else
      raise "Unsupported language!"
    end

    pop_env( :array_start )
    kernel.procedure = p
    return kernel

  end #def BOAST::compute_add_teleseismic_sources_kernel

end
