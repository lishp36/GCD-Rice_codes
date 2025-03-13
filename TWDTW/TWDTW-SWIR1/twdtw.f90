! Fortran 版 TWDTW 算法。
! 编译指令 gfortran ./twdtw.f90 -Ofast -shared -fPIC -o ./twdtw.so
! 根据输入数据类型的不同，分为四个版本。
! double64：时间序列为 real(8)，日期序列为 interger(8)
! double32：时间序列为 real(8)，日期序列为 interger(4)
! single64：时间序列为 real(4)，日期序列为 interger(8)
! single32：时间序列为 real(4)，日期序列为 interger(4)
! 算法来自 R 语言，由董洁改为 MATLAB，由沈若缺改为 Fortran。
! VERSION: v2022.04.14

! TWDTW 核心函数
module dtw_core
   implicit none
   public

   interface dist
      module procedure :: s_dist, d_dist
   end interface dist

   interface multi_dist
      module procedure :: s_multi_dist, d_multi_dist
   end interface multi_dist

   interface multi_tw_dist
      module procedure :: s_multi_tw_dist, d_multi_tw_dist
   end interface multi_tw_dist

   interface tw_dist
      module procedure :: s_tw_dist, d_tw_dist
   end interface tw_dist

   interface tw_gaussian
      module procedure :: s_tw_gaussian, d_tw_gaussian
   end interface tw_gaussian

   interface multi_tw_gaussian
      module procedure :: s_multi_tw_gaussian, d_multi_tw_gaussian
   end interface multi_tw_gaussian

   interface tw_m_dist
      module procedure :: s_tw_m_dist, d_tw_m_dist
   end interface tw_m_dist

   interface tw
      module procedure :: tw_single32, tw_single64, tw_double32, tw_double64
   end interface tw

   interface stw
      module procedure :: stw_single32, stw_single64, stw_double32, stw_double64
   end interface stw

   interface tw_sigma
      module procedure :: s_tw_sigma, d_tw_sigma
   end interface tw_sigma

   interface core
      module procedure :: s_core, d_core
   end interface core

   contains

   subroutine s_dist(x, m, std, n, d)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n)
      real(4), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j))
         end do
      end do
   end subroutine s_dist

   subroutine d_dist(x, m, std, n, d)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n)
      real(8), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j))
         end do
      end do
   end subroutine d_dist

   subroutine s_multi_dist(x, m, t, n, nband, weight, d)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(4), intent(in)  :: x(m, nband), t(n, nband)
      real(4), intent(in)  :: weight(nband)
      real(4), intent(out) :: d(m, n)
      integer :: i, j, iband
      do j = 1, n
         do i = 1, m
            d(i, j) = 0.0_4
            do iband = 1, nband
               d(i, j) = d(i, j) + (abs(x(i, iband) - t(j, iband)) * weight(iband))
            end do
         end do
      end do
   end subroutine s_multi_dist

   subroutine d_multi_dist(x, m, t, n, nband, weight, d)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(8), intent(in)  :: x(m, nband), t(n, nband)
      real(8), intent(in)  :: weight(nband)
      real(8), intent(out) :: d(m, n)
      integer :: i, j, iband
      do j = 1, n
         do i = 1, m
            d(i, j) = 0.0_8
            do iband = 1, nband
               d(i, j) = d(i, j) + (abs(x(i, iband) - t(j, iband)) * weight(iband))
            end do
         end do
      end do
   end subroutine d_multi_dist

   subroutine s_multi_tw_dist(x, m, t, n, nband, weight, tweight, d)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(4), intent(in)  :: x(m, nband), t(n, nband), tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4), intent(out) :: d(m, n)
      integer :: i, j, iband
      do j = 1, n
         do i = 1, m
            d(i, j) = tweight(i, j)
            do iband = 1, nband
               d(i, j) = d(i, j) + (abs(x(i, iband) - t(j, iband)) * weight(iband))
            end do
         end do
      end do
   end subroutine s_multi_tw_dist

   subroutine d_multi_tw_dist(x, m, t, n, nband, weight, tweight, d)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(8), intent(in)  :: x(m, nband), t(n, nband), tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8), intent(out) :: d(m, n)
      integer :: i, j, iband
      do j = 1, n
         do i = 1, m
            d(i, j) = tweight(i, j)
            do iband = 1, nband
               d(i, j) = d(i, j) + (abs(x(i, iband) - t(j, iband)) * weight(iband))
            end do
         end do
      end do
   end subroutine d_multi_tw_dist

   subroutine d_tw_dist(x, m, std, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), tweight(m, n)
      real(8), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) + tweight(i, j)
         end do
      end do
   end subroutine d_tw_dist

   subroutine s_tw_dist(x, m, std, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), tweight(m, n)
      real(4), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) + tweight(i, j)
         end do
      end do
   end subroutine s_tw_dist

   subroutine d_tw_m_dist(x, m, std, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), tweight(m, n)
      real(8), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) * tweight(i, j)
         end do
      end do
   end subroutine d_tw_m_dist

   subroutine s_tw_m_dist(x, m, std, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), tweight(m, n)
      real(4), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) * tweight(i, j)
         end do
      end do
   end subroutine s_tw_m_dist

   subroutine tw_single32(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(4), intent(in)  :: t(m), t_std(n)
      real(4) :: time_diff, alpha, beta
      real(4), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = abs(t(i) - t_std(j))
            tweight(i, j) = 1.0_4 / (1.0_4 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_single32

   subroutine tw_double32(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(4), intent(in)  :: t(m), t_std(n)
      real(8) :: time_diff, alpha, beta
      real(8), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = abs(t(i) - t_std(j))
            tweight(i, j) = 1.0_8 / (1.0_8 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_double32

   subroutine tw_single64(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(4) :: time_diff, alpha, beta
      real(4), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = real(abs(t(i) - t_std(j)), kind=4)
            tweight(i, j) = 1.0_4 / (1.0_4 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_single64

   subroutine tw_double64(t, m, t_std, n, alpha, beta, tweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(8) :: time_diff, alpha, beta
      real(8), intent(out) :: tweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            time_diff = real(abs(t(i) - t_std(j)), kind=8)
            tweight(i, j) = 1.0_8 / (1.0_8 + exp(-alpha * (time_diff - beta)))
         end do
      end do
   end subroutine tw_double64

   subroutine stw_single32(t, k, m, t_std, k_std, n, alpha, beta, stweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(4), intent(in)  :: t(m), t_std(n)
      real(4), intent(in)  :: k(m), k_std(n)
      real(4) :: diff, alpha, beta
      real(4), intent(out) :: stweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            diff = abs(t(i) - t_std(j)) * abs(k(i) - k_std(j))
            stweight(i, j) = 1.0_4 / (1.0_4 + exp(-alpha * (diff - beta)))
         end do
      end do
   end subroutine stw_single32

   subroutine stw_double32(t, k, m, t_std, k_std, n, alpha, beta, stweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(4), intent(in)  :: t(m), t_std(n)
      real(8), intent(in)  :: k(m), k_std(n)
      real(8) :: diff, alpha, beta
      real(8), intent(out) :: stweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            diff = abs(t(i) - t_std(j)) * abs(k(i) - k_std(j))
            stweight(i, j) = 1.0_8 / (1.0_8 + exp(-alpha * (diff - beta)))
         end do
      end do
   end subroutine stw_double32

   subroutine stw_single64(t, k, m, t_std, k_std, n, alpha, beta, stweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(4), intent(in)  :: k(m), k_std(n)
      real(4) :: diff, alpha, beta
      real(4), intent(out) :: stweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            diff = real(abs(t(i) - t_std(j)), kind=4) * abs(k(i) - k_std(j))
            stweight(i, j) = 1.0_4 / (1.0_4 + exp(-alpha * (diff - beta)))
         end do
      end do
   end subroutine stw_single64

   subroutine stw_double64(t, k, m, t_std, k_std, n, alpha, beta, stweight)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(8), intent(in)  :: k(m), k_std(n)
      real(8) :: diff, alpha, beta
      real(8), intent(out) :: stweight(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            diff = abs(t(i) - t_std(j)) *  abs(k(i) - k_std(j))
            stweight(i, j) = 1.0_8 / (1.0_8 + exp(-alpha * (diff - beta)))
         end do
      end do
   end subroutine stw_double64

   subroutine d_tw_gaussian(x, m, std, sigma, n, tweight, d)
      implicit none
      real(8), parameter :: PI = 3.14159265358_8
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(8), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = (sqrt(2 * PI) * sigma(j)) / &
              exp(-(x(i) - std(j)) ** 2 / (2 * (sigma(j) ** 2))) + tweight(i, j)
         end do
      end do
   end subroutine d_tw_gaussian

   subroutine s_tw_gaussian(x, m, std, sigma, n, tweight, d)
      implicit none
      real(4), parameter :: PI = 3.14159265358_4
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(4), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = (sqrt(2 * PI) * sigma(j)) / &
              exp(-(x(i) - std(j)) ** 2 / (2 * (sigma(j) ** 2))) + tweight(i, j)
         end do
      end do
   end subroutine s_tw_gaussian

   subroutine s_multi_tw_gaussian(x, m, t, sigma, n, nband, weight, tweight, d)
      implicit none
      real(4), parameter :: PI = 3.14159265358_4
      integer, intent(in)  :: m, n, nband
      real(4), intent(in)  :: x(m, nband), t(n, nband), sigma(n, nband), tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4), intent(out) :: d(m, n)
      integer :: i, j, iband
      do j = 1, n
         do i = 1, m
            d(i, j) = tweight(i, j)
            do iband = 1, nband
               d(i, j) = d(i, j) + weight(nband) * ((sqrt(2 * PI) * sigma(j, nband)) / &
                 exp(-(x(i, nband) - t(j, nband)) ** 2 / (2 * (sigma(j, nband) ** 2))))
            end do
         end do
      end do
   end subroutine s_multi_tw_gaussian

   subroutine d_multi_tw_gaussian(x, m, t, sigma, n, nband, weight, tweight, d)
      implicit none
      real(8), parameter :: PI = 3.14159265358_8
      integer, intent(in)  :: m, n, nband
      real(8), intent(in)  :: x(m, nband), t(n, nband), sigma(n, nband), tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8), intent(out) :: d(m, n)
      integer :: i, j, iband
      do j = 1, n
         do i = 1, m
            d(i, j) = tweight(i, j)
            do iband = 1, nband
               d(i, j) = d(i, j) + weight(nband) * ((sqrt(2 * PI) * sigma(j, nband)) / &
                 exp(-(x(i, nband) - t(j, nband)) ** 2 / (2 * (sigma(j, nband) ** 2))))
            end do
         end do
      end do
   end subroutine d_multi_tw_gaussian

   subroutine d_tw_sigma(x, m, std, sigma, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(8), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) / sigma(j) / 10 + tweight(i, j)
         end do
      end do
   end subroutine d_tw_sigma

   subroutine s_tw_sigma(x, m, std, sigma, n, tweight, d)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(4), intent(out) :: d(m, n)
      integer :: i, j
      do j = 1, n
         do i = 1, m
            d(i, j) = abs(x(i) - std(j)) / sigma(j) / 10 + tweight(i, j)
         end do
      end do
   end subroutine s_tw_sigma

   real(4) function s_core(d, m, n)
      implicit none
      integer, intent(in)  :: m, n
      real(4) :: d(m, n)
      integer :: i, j
      do i = 2, m
         d(i, 1) = d(i, 1) + d(i-1, 1)
      end do
      do j = 2, n
         d(1, j) = d(1, j) + d(1, j-1)
         do i = 2, m
            d(i, j) = d(i, j) + min(d(i-1, j), d(i-1, j-1), d(i, j-1))
         end do
      end do
      s_core = d(m, n)
   end function s_core

   real(8) function d_core(d, m, n)
      implicit none
      integer, intent(in)  :: m, n
      real(8) :: d(m, n)
      integer :: i, j
      do i = 2, m
         d(i, 1) = d(i, 1) + d(i-1, 1)
      end do
      do j = 2, n
         d(1, j) = d(1, j) + d(1, j-1)
         do i = 2, m
            d(i, j) = d(i, j) + min(d(i-1, j), d(i-1, j-1), d(i, j-1))
         end do
      end do
      d_core = d(m, n)
   end function d_core

end module dtw_core


! Fortran 版 分步骤 TW 与 DTW。
module tw_dtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, std, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), tweight(m, n)
      real(8) :: d(m, n)
      call tw_dist(x, m, std, n, tweight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, std, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), tweight(m, n)
      real(4) :: d(m, n)
      call tw_dist(x, m, std, n, tweight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      double8 = core(d, m, n)
   end function double8

   real(8) function double16(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single8(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      single8 = core(d, m, n)
   end function single8

   real(4) function single16(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_dist(x / factor, m, std, n, tweight, d)
      single16 = core(d, m, n)
   end function single16

end module tw_dtw

! Fortran 版 TWDTW 算法。
module twdtw
   use dtw_core
   use tw_dtw
   implicit none
   public
   contains

   real(8) function double64(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(8), intent(in)  :: x(m), std(n)
      real(8) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      double64 = doublex(x, m, std, n, tweight)
   end function double64

   real(8) function double32(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer, intent(in)  :: t(m), t_std(n)
      real(8), intent(in)  :: x(m), std(n)
      real(8) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      double32 = doublex(x, m, std, n, tweight)
   end function double32


   real(4) function single64(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer(8), intent(in)  :: t(m), t_std(n)
      real(4), intent(in)  :: x(m), std(n)
      real(4) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      single64 = singlex(x, m, std, n, tweight)
   end function single64

   real(4) function single32(x, t, m, std, t_std, n, alpha, beta)
      implicit none
      integer, intent(in)  :: m, n
      integer, intent(in)  :: t(m), t_std(n)
      real(4), intent(in)  :: x(m), std(n)
      real(4) :: tweight(m, n), alpha, beta
      call tw(t, m, t_std, n, alpha, beta, tweight)
      single32 = singlex(x, m, std, n, tweight)
   end function single32

end module twdtw

! Fortran 版 DTW 算法。
module dtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, std, n)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n)
      real(8) :: d(m, n)
      call dist(x, m, std, n, d)
      doublex = core(d, m, n)
   end function doublex


   real(4) function singlex(x, m, std, n)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n)
      real(4) :: d(m, n)
      call dist(x, m, std, n, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in) :: x(m)
      real(8), intent(in)  :: std(n), factor
      real(8) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      double8 = core(d, m, n)
   end function double8

   real(4) function single8(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in) :: x(m)
      real(4), intent(in)  :: std(n), factor
      real(4) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      single8 = core(d, m, n)
   end function single8

   real(8) function double16(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in) :: x(m)
      real(8), intent(in)  :: std(n), factor
      real(8) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single16(x, m, std, n, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in) :: x(m)
      real(4), intent(in)  :: std(n), factor
      real(4) :: d(m, n)
      call dist(x / factor, m, std, n, d)
      single16 = core(d, m, n)
   end function single16

end module dtw

! Fortran 版 多元 DTW 算法。
module multidtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, t, n, nband, weight)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(8), intent(in)  :: x(m, nband), t(n, nband)
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_dist(x, m, t, n, nband, weight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, t, n, nband, weight)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(4), intent(in)  :: x(m, nband), t(n, nband)
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_dist(x, m, t, n, nband, weight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, t, n, nband, weight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(1), intent(in) :: x(m, nband)
      real(8), intent(in)  :: t(n, nband), factor
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_dist(x / factor, m, t, n, nband, weight, d)
      double8 = core(d, m, n)
   end function double8

   real(4) function single8(x, m, t, n, nband, weight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(1), intent(in) :: x(m, nband)
      real(4), intent(in)  :: t(n, nband), factor
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_dist(x / factor, m, t, n, nband, weight, d)
      single8 = core(d, m, n)
   end function single8

   real(8) function double16(x, m, t, n, nband, weight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(2), intent(in) :: x(m, nband)
      real(8), intent(in)  :: t(n, nband), factor
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_dist(x / factor, m, t, n, nband, weight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single16(x, m, t, n, nband, weight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(2), intent(in) :: x(m, nband)
      real(4), intent(in)  :: t(n, nband), factor
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_dist(x / factor, m, t, n, nband, weight, d)
      single16 = core(d, m, n)
   end function single16

end module multidtw

! Fortran 版 多元 TWDTW 算法，步骤 TW 与 DTW。
module multi_tw_dtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, t, n, nband, weight, tweight)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(8), intent(in)  :: x(m, nband), t(n, nband), tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_tw_dist(x, m, t, n, nband, weight, tweight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, t, n, nband, weight, tweight)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(4), intent(in)  :: x(m, nband), t(n, nband), tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_tw_dist(x, m, t, n, nband, weight, tweight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, t, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(1), intent(in) :: x(m, nband)
      real(8), intent(in)  :: t(n, nband), factor, tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_tw_dist(x / factor, m, t, n, nband, weight, tweight, d)
      double8 = core(d, m, n)
   end function double8

   real(4) function single8(x, m, t, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(1), intent(in) :: x(m, nband)
      real(4), intent(in)  :: t(n, nband), factor, tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_tw_dist(x / factor, m, t, n, nband, weight, tweight, d)
      single8 = core(d, m, n)
   end function single8

   real(8) function double16(x, m, t, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(2), intent(in) :: x(m, nband)
      real(8), intent(in)  :: t(n, nband), factor, tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_tw_dist(x / factor, m, t, n, nband, weight, tweight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single16(x, m, t, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(2), intent(in) :: x(m, nband)
      real(4), intent(in)  :: t(n, nband), factor, tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_tw_dist(x / factor, m, t, n, nband, weight, tweight, d)
      single16 = core(d, m, n)
   end function single16

end module multi_tw_dtw

! Fortran 版 分步骤 TW 与 DTW，中间用乘法。
module tw_m_dtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, std, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), tweight(m, n)
      real(8) :: d(m, n)
      call tw_m_dist(x, m, std, n, tweight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, std, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), tweight(m, n)
      real(4) :: d(m, n)
      call tw_m_dist(x, m, std, n, tweight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_m_dist(x / factor, m, std, n, tweight, d)
      double8 = core(d, m, n)
   end function double8

   real(8) function double16(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_m_dist(x / factor, m, std, n, tweight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single8(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_m_dist(x / factor, m, std, n, tweight, d)
      single8 = core(d, m, n)
   end function single8

   real(4) function single16(x, m, std, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_m_dist(x / factor, m, std, n, tweight, d)
      single16 = core(d, m, n)
   end function single16

end module tw_m_dtw

! Fortran 版 分步骤 TW 与 DTW，改用正态分布的概率密度。
module tw_dtw_gaussian
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, std, sigma, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(8) :: d(m, n)
      call tw_gaussian(x, m, std, sigma, n, tweight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, std, sigma, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(4) :: d(m, n)
      call tw_gaussian(x, m, std, sigma, n, tweight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_gaussian(x / factor, m, std, sigma, n, tweight, d)
      double8 = core(d, m, n)
   end function double8

   real(8) function double16(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_gaussian(x / factor, m, std, sigma, n, tweight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single8(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_gaussian(x / factor, m, std, sigma, n, tweight, d)
      single8 = core(d, m, n)
   end function single8

   real(4) function single16(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_gaussian(x / factor, m, std, sigma, n, tweight, d)
      single16 = core(d, m, n)
   end function single16

end module tw_dtw_gaussian

! Fortran 版 多元 TWDTW 算法，步骤 TW 与 DTW，改用正态分布的概率密度。
module multi_tw_dtw_gaussian
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, t, sigma, n, nband, weight, tweight)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(8), intent(in)  :: x(m, nband), t(n, nband), sigma(n, nband), tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_tw_gaussian(x, m, t, sigma, n, nband, weight, tweight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, t, sigma, n, nband, weight, tweight)
      implicit none
      integer, intent(in)  :: m, n, nband
      real(4), intent(in)  :: x(m, nband), t(n, nband), sigma(n, nband), tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_tw_gaussian(x, m, t, sigma, n, nband, weight, tweight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, t, sigma, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(1), intent(in) :: x(m, nband)
      real(8), intent(in)  :: t(n, nband), sigma(n, nband), factor, tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_tw_gaussian(x / factor, m, t, sigma, n, nband, weight, tweight, d)
      double8 = core(d, m, n)
   end function double8

   real(4) function single8(x, m, t, sigma, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(1), intent(in) :: x(m, nband)
      real(4), intent(in)  :: t(n, nband), sigma(n, nband), factor, tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_tw_gaussian(x / factor, m, t, sigma, n, nband, weight, tweight, d)
      single8 = core(d, m, n)
   end function single8

   real(8) function double16(x, m, t, sigma, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(2), intent(in) :: x(m, nband)
      real(8), intent(in)  :: t(n, nband), sigma(n, nband), factor, tweight(m, n)
      real(8), intent(in)  :: weight(nband)
      real(8) :: d(m, n)
      call multi_tw_gaussian(x / factor, m, t, sigma, n, nband, weight, tweight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single16(x, m, t, sigma, n, nband, weight, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n, nband
      integer(2), intent(in) :: x(m, nband)
      real(4), intent(in)  :: t(n, nband), sigma(n, nband), factor, tweight(m, n)
      real(4), intent(in)  :: weight(nband)
      real(4) :: d(m, n)
      call multi_tw_gaussian(x / factor, m, t, sigma, n, nband, weight, tweight, d)
      single16 = core(d, m, n)
   end function single16

end module multi_tw_dtw_gaussian

! Fortran 版 分步骤 TW 与 DTW，距离除以标准差。
module tw_dtw_sigma
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, m, std, sigma, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(8), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(8) :: d(m, n)
      call tw_sigma(x, m, std, sigma, n, tweight, d)
      doublex = core(d, m, n)
   end function doublex

   real(4) function singlex(x, m, std, sigma, n, tweight)
      implicit none
      integer, intent(in)  :: m, n
      real(4), intent(in)  :: x(m), std(n), sigma(n), tweight(m, n)
      real(4) :: d(m, n)
      call tw_sigma(x, m, std, sigma, n, tweight, d)
      singlex = core(d, m, n)
   end function singlex

   real(8) function double8(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_sigma(x / factor, m, std, sigma, n, tweight, d)
      double8 = core(d, m, n)
   end function double8

   real(8) function double16(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(8), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(8) :: d(m, n)
      call tw_sigma(x / factor, m, std, sigma, n, tweight, d)
      double16 = core(d, m, n)
   end function double16

   real(4) function single8(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(1), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_sigma(x / factor, m, std, sigma, n, tweight, d)
      single8 = core(d, m, n)
   end function single8

   real(4) function single16(x, m, std, sigma, n, tweight, factor)
      implicit none
      integer, intent(in)  :: m, n
      integer(2), intent(in)  :: x(m)
      real(4), intent(in)  :: std(n), sigma(n), tweight(m, n), factor
      real(4) :: d(m, n)
      call tw_sigma(x / factor, m, std, sigma, n, tweight, d)
      single16 = core(d, m, n)
   end function single16

end module tw_dtw_sigma

! STWDTW (Shape and Time Weighted Dynamic Time Warping)
module multi_stwdtw
   use dtw_core
   implicit none
   public
   contains

   real(8) function doublex(x, t, k, m, std, t_std, k_std, n, nband, weight, alpha, beta)
      implicit none
      integer, intent(in) :: m, n, nband, t(m), t_std(n)
      real(8), intent(in) :: x(m, nband), k(m, nband)
      real(8), intent(in) :: std(n, nband), k_std(n, nband)
      real(8), intent(in) :: weight(nband), alpha, beta
      real(8) :: d(m, n), stweight(m, n), d_all(m, n)
      integer :: iband
      d_all = 0
      do iband = 1, nband
         call stw(t, k(:, iband), m, t_std, k_std(:, iband), n, alpha, beta, stweight)
         call dist(x(:, iband), m, std(:, iband), n, d)
         d_all = d_all + d * stweight * weight(iband)
      end do
      doublex = core(d_all, m, n)
   end function doublex

   real(4) function singlex(x, t, k, m, std, t_std, k_std, n, nband, weight, alpha, beta)
      implicit none
      integer, intent(in) :: m, n, nband, t(m), t_std(n)
      real(4), intent(in) :: x(m, nband), k(m, nband)
      real(4), intent(in) :: std(n, nband), k_std(n, nband)
      real(4), intent(in) :: weight(nband), alpha, beta
      real(4) :: d(m, n), stweight(m, n), d_all(m, n)
      integer :: iband
      d_all = 0
      do iband = 1, nband
         call stw(t, k(:, iband), m, t_std, k_std(:, iband), n, alpha, beta, stweight)
         call dist(x(:, iband), m, std(:, iband), n, d)
         d_all = d_all + d * stweight * weight(iband)
      end do
      singlex = core(d_all, m, n)
   end function singlex

end module multi_stwdtw
