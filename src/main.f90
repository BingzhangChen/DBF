PROGRAM MAIN
USE MPI
IMPLICIT NONE

real,    parameter :: pi      = 3.1415926535
integer, parameter :: d_per_y = 360   !How many days a year
integer, parameter :: s_per_d = 86400 !How many seconds a day

!Number of days running:
integer, parameter :: NDAYS   = d_per_y*5 !5 years

!Timestep (unit: day)
real, parameter    :: dt = 300.d0/dble(s_per_d)

!Number of integration:
integer :: NI=100

!Tracers contain NO3 and PHY
!Initialize tracers:
real :: NO3 = 1., PHY=1D-3
real :: NO3_avg = 0.  !Average NO3 for saving data
real, parameter :: mort   = .1d0  !Mortality rate for phyto. unit: d-1 (mmolN L-1)-1

real :: Pmort=0.  !Phytoplankton mortality for a single species
real :: TPmort=0.  !Total Phytoplankton mortality loss that goes back to NO3

!Amplitude of nutrient inflow
integer, parameter :: NAm = 9
integer, parameter :: Nfreq = 9
real, parameter :: Am(NAm) = [0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0,0.8d0, 0.9d0] 
real, parameter :: Freq(Nfreq) = [1.d0, 2d0, 4d0, 8d0, 16d0, 32d0, 64d0, 128d0, 256d0] 

!Nutrient turnover rate
real, parameter :: Di = .2d0 !Unit: per day
real            :: Di_ = 0d0 !Time-dependent flushing rate

!Inflow constant of nutrient:
real, parameter :: N0=1.d0 !Nutrient concentration of nutrient reservoire

!diversity controlling parameter
real, parameter :: phi = .5d0 !Larger phi leads to greater diversity

real    :: cur_day ! Current simulation time (day)
integer :: i,j,k,jj,v,kk

!Species richness from 2 to 10
integer, parameter :: M = 10, NR=10 !Number of repetitions

integer, parameter :: TN= 100  ! Total number of size classes
real :: LNV(TN)
real, allocatable  :: LNV_(:), PHY_(:), PHY_avg(:), NPP_(:), NPP_avg(:)

!File names for monoculture simulations
character(Len=30), allocatable  :: mfname(:)

!Vectors  for saving NO3, PHY, and PP for each simulation of monocultures
real, allocatable :: mNO3(:), mPHY(:), mPP(:)
real, allocatable :: mNO3_avg(:), mPHY_avg(:), mPP_avg(:)
integer :: AllocateStatus = 0

!Growth rate at the mean size, first and second derivatives of growth
real :: mu, dx
real, external :: logV

real, parameter :: Vmin = -2.726 !0.5 micron
real, parameter :: Vmax= 13.168 !100 micron

!Variables for timing
REAL         :: T1, T2

!Filename
character(LEN=20) :: fname = 'Mod.out'

!Scratch character
character(LEN=9) :: str = 'MOD'

!Counter for averaging
integer :: count = 0

! MPI variables:
integer :: ierr, ntasks

!End of declaration

!Start 
CALL CPU_TIME(T1)

ntasks=NR

! ***** Initialize MPI *****
call MPI_INIT(ierr)

! Returns the size of the group associated with a communicator.
call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

! Determines the rank of the calling process in the communicator
call MPI_COMM_RANK(MPI_COMM_WORLD, k, ierr)

! Notes: each mpi process simulates a random replication

! Number of integrations:
NI=NDAYS*INT(1.d0/dt)

call random_seed()

!Random sampling for 10 combinations
DO i = 2, M !Iteration of richness
      !Generate a random number between 1 and TN
      allocate(LNV_(i), stat = allocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating LNV_ ***"
      LNV_(:) = 0.d0

      do j = 1, i   ! Totally i species, same size for different disturbance (Am)
         call random_number(dx)
         LNV_(j) = dx*(Vmax - Vmin) + Vmin  !Randomly select i species
      enddo
 
      allocate(PHY_(i),  stat = allocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating PHY_ ***"
      PHY_(:) = 0.d0
     
      allocate(PHY_avg(i),  stat = allocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating PHY_avg ***"
      PHY_avg(:) = 0.d0

      allocate(NPP_(i),  stat = allocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating NPP_ ***"
      NPP_(:) = 0.d0

      allocate(NPP_avg(i),  stat = allocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating NPP_avg ***"
      NPP_avg(:) = 0.d0

      !Initialise variables for monocultures:
      allocate(mNO3(i), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating mNO3 ***"
      mNO3(:)=NO3

      allocate(mNO3_avg(i), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating mNO3_avg ***"
      mNO3_avg = mNO3

      allocate(mPHY(i), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating mPHY ***"
      mPHY(:)=PHY

      allocate(mPHY_avg(i), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating mPHY_avg ***"
      mPHY_avg=mPHY

      allocate(mPP(i), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating mPP ***"
      mPP(:)=0d0

      allocate(mPP_avg(i), stat=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Problem in allocating mPP_avg ***"
      mPP_avg=mPP

      DO v = 1, NAm !run the same species under different Am
        DO kk = 1, Nfreq

          !Initialise variables for polycultures:
          NO3 = 1.d0; PHY=1D-3
          NO3_avg = NO3
     
          do j = 1, i   ! Totally i species, same initial conditions for different disturbance (Am)
             PHY_(j) = PHY/dble(i)
          enddo
     
          PHY_avg(:) = PHY_(:)

          !Calculate initial growth rate
          do jj = 1,i 
             call growth(LNV_(jj), NO3, mu) 
             NPP_(jj) = PHY_(jj)*mu
          enddo

          NPP_avg = NPP_

          !Create file name for polycultures:
          write(fname,100) 'P',v,'_', i,'_',k,'_',kk,'.out'
     
          !saving results to polyculture file
          open (unit=10, file = fname, status = 'replace')
          write(10, 101) 'Day', 'NO3', LNV_, LNV_
     
          !Record initial conditions:
          write(10, 102) 0., NO3, PHY_, NPP_avg
          close(10)

          !Create file name for monocultures:
          allocate(mfname(i),  stat=AllocateStatus)
          IF (AllocateStatus /= 0) STOP "*** Problem in allocating mfname ***"
          mfname(:) = 'M.out'

         do j = 1, i
           write(mfname(j),103) 'M',v,'_', i,'_',k,'_',kk,'_',j,'.out'

           !Writing headers to monoculture files
           open (unit=10, file = mfname(j), status = 'replace')
           write(10, 104) 'Day', 'NO3', LNV_(j), 'NPP'
     
           !Record initial conditions:
           mNO3(:)=NO3
           mPHY(:)=PHY

           !Calculate initial growth rate
           do jj = 1,i 
              call growth(LNV_(jj), NO3, mu) 
              mPP(jj) = PHY_(jj)*mu
           enddo
           mPHY_avg = mPHY
           mPP_avg = mPP

           write(10, 102) 0., mNO3(j), mPHY(j), mPP(j)
           close(10)
         enddo
    
         Do j = 1, NI  !Start integration

            !Update counter
            count = count + 1

            !Current time (in days)
            cur_day = dt*dble(j)

            !Calculate flushing rate
            Di_= Di*(1.d0+Am(v)*sin(2.d0*pi*cur_day*Freq(kk)/dble(d_per_y))) 
   
            !Save data every day:
            If (mod(j,int(1./dt)) .eq. 0) then

              !Calculate mean values of nutrient and phyto. biomass during this period
              !Write polyculture data into the file

              !Compute average values throughout this period
              open (unit=10, file = fname, status = 'old', action='write', position='append')
              write(10, 102) cur_day, NO3_avg, PHY_avg, NPP_avg
              close(10)

              !Reset average values of polycultures
              NO3_avg = 0d0
              NPP_avg(:) = 0d0
              PHY_avg(:) = 0d0

              !Monocultures

              !Write monoculture data into the file
              do jj = 1, i
                 open (unit=10, file = mfname(jj), status = 'old', action='write', position='append')
                 write(10, 102) cur_day, mNO3_avg(jj), mPHY_avg(jj), mPP_avg(jj)
                 close(10)
              enddo

              !Reset average values of monocultures
              mNO3_avg(:) = 0d0
              mPP_avg(:) = 0d0
              mPHY_avg(:) = 0d0

              !Reset counter
              count = 0
            Endif !End of saving data to external files
     
            !Calculate growth rate for each species (polyculture simulations)
            PHY = sum(PHY_(:))
            TPmort=0d0
            do jj = 1,i 
               call growth(LNV_(jj), NO3, mu) 

               !Production of each species
               NPP_(jj)= PHY_(jj)*mu
               Pmort = mort*PHY**(1d0-phi)*PHY_(jj)**(1d0+phi) !Use the approach of Record et al. (2014) to sustain phyto. diversity
               PHY_(jj) = PHY_(jj) + ( (mu- Di_)*PHY_(jj)-Pmort )*dt
               PHY_(jj) = max(PHY_(jj), 0d0) !Ensure positive number
               TPmort   = TPmort + Pmort

               !Update PHY_avg
               PHY_avg(jj) = (PHY_avg(jj)*dble(count) + PHY_(jj))/dble(count+1)

               !Update PP_avg
               NPP_avg(jj) = (NPP_avg(jj)*dble(count) + NPP_(jj))/dble(count+1)
            enddo
  
            !Iteration of NO3 (Phytoplankton mortality goes back to nutrient)
            NO3 = NO3 + dt*( (N0-NO3)*Di_ - sum(NPP_(:)) + TPmort)

            !Update NO3_avg
            NO3_avg = (NO3_avg*dble(count) + NO3)/dble(count+1)

            !Simulate monocultures
            do jj = 1,i 
               call growth(LNV_(jj), mNO3(jj), mu) 
               mPP(jj)= mPHY(jj)*mu
               Pmort   = mort*mPHY(jj)**2
               mPHY(jj) = mPHY(jj) + ((mu - Di_)*mPHY(jj) - Pmort)*dt
               mNO3(jj) = mNO3(jj)*(1.d0+ (N0 - mNO3(jj))*Di*dt) + (- mPP(jj) + Pmort)*dt

               !Update mNO3_avg
               mNO3_avg(jj) = (mNO3_avg(jj)*dble(count) + mNO3(jj))/dble(count+1)

               !Update mPHY_avg
               mPHY_avg(jj) = (mPHY_avg(jj)*dble(count) + mPHY(jj))/dble(count+1)

               !Update mPP_avg
               mPP_avg(jj) = (mPP_avg(jj)*dble(count) + mPP(jj))/dble(count+1)

            enddo

          Enddo !End of timestep integration
        DEALLOCATE(mfname)
      ENDDO !End of loop of frequencies
     ENDDO  !End of loop of Ams
     DEALLOCATE(LNV_)
     DEALLOCATE(PHY_)
     DEALLOCATE(PHY_avg)
     DEALLOCATE(NPP_)
     DEALLOCATE(NPP_avg)
     DEALLOCATE(mNO3)
     DEALLOCATE(mNO3_avg)
     DEALLOCATE(mPP)
     DEALLOCATE(mPP_avg)
     DEALLOCATE(mPHY)
     DEALLOCATE(mPHY_avg)
ENDDO !End of richness gradients

! Synchronize all processes:
call MPI_BARRIER (MPI_COMM_WORLD,ierr)

!End mpi
call MPI_finalize(ierr)

CALL CPU_TIME(T2)
WRITE(6, '("TIME = ",F8.3," hours.")') (T2-T1)/3600.0

100 format(4(A1,I0),A4) 
101 format(6x,2(A3, 6x), 20(1pe12.4, 1x))
102 format(F12.3, 6x,  60(1pe12.2, 1x))   
103 format(5(A1,I0),A4) 
104 format(6x,2(A3, 6x), 1pe12.4, 3x, A3)
END PROGRAM MAIN

subroutine growth(PMU, N,mu) 
implicit none
real, intent(in)  :: PMU, N
real, intent(out) :: mu
real, parameter   :: mu0 = 1.d0, alphamu = .2, betamu = 0.d0
real, parameter   :: K0N = .5, alphaK = .3
real :: mu0hat
real :: Kn, fN

    mu0hat = mu0*exp(alphamu*PMU + betamu*PMU**2)
    Kn = K0N*exp(alphaK*PMU)
    fN = N/(Kn+N)
    mu = mu0hat*fN
end subroutine growth

pure real function logV(ESD)
implicit none
real, intent(in) :: ESD
real, parameter  :: pi = 3.14159265357989
logV = log(pi/6.*ESD**3)
end function logV

!Function calculating biomass-weighted mean trait and trait variance
subroutine Mean_VAR(N, L_,P_, avgL, VarL) 
implicit none
integer, intent(in) :: N !Number of species
real, intent(in) :: L_(N) !Vector of trait values
real, intent(in) :: P_(N) !Vector of species biomass
real, intent(out):: avgL !Mean Trait
real, intent(out):: VarL !Trait variance
integer :: j

avgL = 0.d0
VarL = 0.d0
do j = 1, N
    avgL = avgL+P_(j)*L_(j)
enddo
avgL = avgL/sum(P_(:))
do j = 1, N
   VarL = VarL + P_(j)*(L_(j) - avgL)**2
enddo
VarL = VarL/sum(P_(:))
return
end subroutine Mean_VAR