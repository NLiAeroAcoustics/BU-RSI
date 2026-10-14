      function FAC(n)

      implicit real*8 (a-h, o-z)
 
      if (n .gt. 1) got o 1
       FAC = 1.D0
        return

 1    FAC = 1.D0
       do 2 j = 1, N
        FAC = FAC * DFLOAT(j)
 2    continue

      return
       end
