      function isig (x)

      implicit real*8 (a-h, o-z)
      integer isig
      
      isig = 1
       if (x .lt. 0.D0) isig = -1
       
      return
       end
