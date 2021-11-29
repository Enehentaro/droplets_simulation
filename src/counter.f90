module counter
	implicit none
	contains

	function count_txt() result(counter)
		character(50),parameter :: f_name='name.txt' !読み込むtxtファイルの名前を付ける
		integer unit
		integer counter
		counter=0

		open(newunit=unit,file=f_name,status='old')
		do
			read(unit,*,end=999)
			counter=counter+1
		end do
		999 close(unit)
	end function count_txt

	subroutine read_txt(num_file,field_name)
		character(50),parameter :: f_name='name.txt' !読み込むtxtファイルの名前を付ける
		integer unit,i
		integer,intent(in) :: num_file
		character(50), allocatable :: field_name(:)
		allocate(field_name(num_file))

		open(newunit=unit,file=f_name,status='old')
			read(unit,*)(field_name(i),i=1,num_file)
		close(unit)
	end subroutine read_txt

end module counter