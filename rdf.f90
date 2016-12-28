		PROGRAM RDF

! CALCULATE RDF										 
! Jun Taek Lee <leejunta@grinnell.edu>

!*******************************USE MODULE************************************************

		use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
	    use xtc

!*******************************ASSIGN VARIABLES******************************************

		implicit none

		real(8), allocatable			::	com(:,:,:)										!Center of mass (molecule A)
		real(8), allocatable			::	car(:,:,:,:)										!Atom Cartesian coordinates	
		real(8), allocatable			::	comb(:,:,:)										!Center of mass (molecule B)
		real(8)							::	l(3), rho1, rho									!length of box, rho for g1, rho for g2
		real(8)							::	deltab											!length of bin
		real(8), allocatable			::	massa(:,:)										!mass of atom
		real(8), allocatable			::	massm(:)										!mass of molecule
		real(8)							::	rij(3), r, rcomp								!displacement between molecules
		real(8)							::	rlower, rupper, ravg							!radius from molecule A
		real(8), allocatable			::	sumg(:), sumg1(:)
		real(8), allocatable			::	histavg1(:), histavg(:), g(:), g1(:)		
		double precision, PARAMETER 	:: 	pi=3.1415926535d+0					
		integer, allocatable			::	nm(:)											!number of molecules
		integer, allocatable			::	na(:)											!number of atoms
		integer, allocatable			::	hist(:,:)	
		integer							::	bin, maxbin				
		integer							::	tstart
		integer							::	cube
		integer							::	maxmolsp, imolsp								!maximum # of molecular species, current	
		integer							::	omol											!arbitrary molecule number
		integer							::	x,y,z											!dimensions x,y,z
		integer							::	tnm, tmaxstep									!total # of molecules, time steps
		integer 						::	iatom, itatom, xyz, itime						!current atom, current atom, dimensions, current time
		integer							::	imol, imola, imolb								!molecule number
		character*40					::	outfile				

 		character*40					:: infile								
   		real, allocatable				:: pos(:,:)
    	integer 						:: NATOMS, STEP, STAT
    	real 							:: box(3,3), prec, time, box_trans(3,3)
    	type(C_PTR) 					:: xd_c
    	type(xdrfile), pointer 			:: xd
    	logical 						:: ex
	
!****************************READ SHELL SCRIPT********************************************

		read(5,*)	infile, outfile, maxmolsp	
		read(5,*)	deltab, tstart, maxbin
!		infile = "ether-g.xtc"//C_NULL_CHAR
		
!****************************OPEN XTC FILE THROUGH MODULE*********************************		
		
    	inquire(file=trim(infile),exist=ex)

    	STAT = read_xtc_natoms(infile,NATOMS)
    	allocate(pos(3,NATOMS))

    	xd_c = xdrfile_open(infile,"r")
    	call c_f_pointer(xd_c,xd)		
    	
!****************************ALLOCATE, READ, CALCULATE VARIABLES**************************

		allocate(na(maxmolsp))									
		allocate(nm(maxmolsp))									
		allocate(massa(maxmolsp,na(maxmolsp)))									
		allocate(massm(maxmolsp))	
	
		massm = 0								
		do imolsp = 1, maxmolsp								
		  read(5,*) nm(imolsp), na(imolsp)		
		end do  		
			
		read(5,*) tnm

		do imolsp = 1, maxmolsp		
		  do itatom = 1,na(imolsp)						
		    read(5,'(F7.3)') massa(imolsp,itatom)	
		  end do
		end do  		

		write(*,*) "Done Reading Mass"

		itime=0
		cube=9
		
		allocate(com(3,maxmolsp,nm(maxmolsp)))
		allocate(comb(3,maxmolsp,cube*nm(maxmolsp)))
		allocate(car(3,maxmolsp,nm(maxmolsp),na(maxmolsp)))
		allocate(hist(maxmolsp,maxbin))
		allocate(histavg(maxbin))
		allocate(histavg1(maxbin))
		allocate(g(maxbin))
		allocate(g1(maxbin))	
		allocate(sumg(maxbin))
		allocate(sumg1(maxbin))
				
!******************************OPEN OUTPUT FILE*******************************************
		
		open(unit=11,file=outfile)								
		
!**************************READ ATOM POSITIONS********************************************	
	
		write(*,*) "Reading xtc"
		write(*,*) "..."
		write(*,*) "..."
		write(*,*) "..."

    	STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
        box = transpose(box_trans)
		do while ( STAT == 0 )
			itatom=0
			if (time < tstart) goto 91
			if(time == tstart) then
				l(1) = 10.0d0*box(1,1)
				l(2) = 10.0d0*box(2,2)
				l(3) = 10.0d0*box(3,3)
				rho1 = nm(1)/(l(1)*l(2)*l(3))
				rho	 = tnm/(l(1)*l(2)*l(3))
			end if
			
!******************************ASSEMBLE SEPARATED MOLECULES TOGETHER (PBC)****************
			do imolsp = 1, maxmolsp
				do imol = 1, nm(imolsp)
					do iatom = 1, na(imolsp)
						itatom = itatom + 1
						car(:,imolsp,imol,iatom)=pos(:,itatom)
						do xyz = 1,3
							if ((car(xyz,imolsp,imol,iatom)-car(xyz,imolsp,imol,1)) > (l(xyz)/20)) then
								car(xyz,imolsp,imol,iatom) = car(xyz,imolsp,imol,iatom) - (l(xyz)/10)
							else if ((car(xyz,imolsp,imol,iatom)-car(xyz,imolsp,imol,1)) < (-l(xyz)/20)) then
								car(xyz,imolsp,imol,iatom) = car(xyz,imolsp,imol,iatom) + (l(xyz)/10)
							end if
						end do
						
!******************************CALCULATE COM FOR EACH MOLECULE****************************

					    com(:,imolsp,imol) = com(:,imolsp,imol) + 10.0d0*car(:,imolsp,imol,iatom)*massa(imolsp,iatom)
						if ((time.eq.tstart) .and. (imol.eq.1)) then	
							massm(imolsp) = massm(imolsp) + massa(imolsp,iatom)
						end if				
					end do							
			    	com(:,imolsp,imol) = com(:,imolsp,imol) / massm(imolsp)	
				end do
			end do
			
!**************************PBC: CREATE 8 MORE BOXES***************************************

			do imolsp=1,maxmolsp
				omol=0
				do imol=1,nm(imolsp)
					do xyz=1,3	
				    	if(com(xyz,imolsp,imol) >= l(xyz)) then
							com(xyz,imolsp,imol) = com(xyz,imolsp,imol) - l(xyz)
						else if(com(xyz,imolsp,imol) < 0) then
							com(xyz,imolsp,imol) = com(xyz,imolsp,imol) + l(xyz)
						end if
					end do
					do x=0,2
						do y=0,2
							do z=0,2
								omol=omol+1
								comb(1,imolsp,omol)=com(1,imolsp,imol)+ (x-1)*l(1)
								comb(2,imolsp,omol)=com(2,imolsp,imol)+ (y-1)*l(2)
								comb(3,imolsp,omol)=com(3,imolsp,imol)+ (z-1)*l(3)
							end do
						end do
					end do
				end do
			end do

!*******************CALCULATE DISTANCE AND ADD TO BINS************************************

			do imolsp=1,maxmolsp
				do imola = 1, nm(1)
					do imolb = 1, (cube*nm(imolsp))
						r=0
						rij(:) = com(:,1,imola) - comb(:,imolsp,imolb)
						do xyz=1,3
							rcomp = rij(xyz)**2
							r = r + rcomp
						end do
						r = sqrt(r)
						if (r==0) cycle
						bin  = floor(r/deltab)+1
						if (bin <= maxbin .and. bin > 1) then
							hist(imolsp,bin) = hist(imolsp,bin) + 1
						end if		
					end do
				end do	
			end do
			itime=itime+1
91			write(*,*) "Check", time
			com=0
    		STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
		end do
		STAT = xdrfile_close(xd)
		write(*,*) "Finished reading"
				
!************************NORMALIZATION****************************************************

		write(*,*) "Calculating"
		write(*,*) "..."
		write(*,*) "..."
		
		sumg=0
		sumg1=0
		do bin=1, maxbin
			rlower					= real(bin-1)*deltab
			rupper					= rlower + deltab
			ravg					= rlower + deltab/2
			histavg1(bin)			= real(hist(1,bin))/(nm(1)*itime)	
			g1(bin)					= histavg1(bin)/(4*pi*ravg**2*deltab*rho1)
			histavg(bin)			= real(hist(1,bin)+hist(2,bin))/(nm(1)*itime)	
			g(bin)					= histavg(bin)/(4*pi*ravg**2*deltab*rho1)
			if (bin>1) then
				sumg1(bin)			= sumg1(bin-1) + g1(bin)
				sumg(bin)			= sumg(bin-1) + g(bin)
			end if
			write(11,'(2x,f8.3,a,2x,f8.3,a,2x,f8.3,a,2x,f8.3,a,2x,f8.3)') ravg,",", g1(bin),",", g(bin),",",sumg1(bin),",",sumg(bin)
		end do

		write(*,*) "Done!"

!******************************DEALLOCATE AND CLOSE OUTPUT FILE***************************
		
		close(unit=11)											
		
		deallocate(na)				
		deallocate(nm)							
		deallocate(massa)										
		deallocate(massm)	
		deallocate(com)
		deallocate(comb)
		deallocate(car)
		deallocate(hist)
		deallocate(histavg)
		deallocate(histavg1)
		deallocate(g)
		deallocate(g1)	
	    deallocate(pos)	
	    deallocate(sumg)
	    deallocate(sumg1)
		
!*****************************END PROGRAM*************************************************
				
		END PROGRAM											
				
				
				
				
				
				
				