!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Nonequilibrium structure: module containing nonequilibrium forces available for MD
!> \author James Reid [JR] - University of Queensland - 06.2017 - Based on thermostat_type
! *****************************************************************************
MODULE nonequilibrium_types
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_constants,                 ONLY: do_ne_colour,&
                                              do_ne_field
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_list_types,             ONLY: particle_list_type
   USE simpar_types,                    ONLY: simpar_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: NE_type, &
             !            allocate_NE,&
             !            release_NE,&
             !            retain_NE,&
             create_NE_type
!            release_NE_type,&
!            NE_info_type,&
!            release_NE_info,&
!            set_NE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'Non-equilibrium_methods'

   ! Single NE_type JR - should trim down unused parts!!
! *****************************************************************************
   TYPE NE_type
      INTEGER :: id_nr, ref_count, type_of_NE
      CHARACTER(LEN=default_string_length) :: label
      REAL(KIND=dp)                            :: FX, FY, FZ
   END TYPE NE_type

CONTAINS

! *****************************************************************************
!> \brief Create a nonequilibrium type (based off of thermostat_type)
!> \param NEtype ...
!> \param simpar ...
!> \param md_section ...
!> \param force_env ...
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino Modified by JR
! **************************************************************************************************
   SUBROUTINE create_NE_type(NEtype, simpar, md_section, force_env)
      TYPE(NE_type), POINTER                             :: NEtype
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(section_vals_type), POINTER                   :: md_section
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'create_NE_type', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ipart, stat
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: explicit, shell_present
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(section_vals_type), POINTER                   :: section

      simpar%do_NEforce = .FALSE.

      section => section_vals_get_subs_vals(md_section, "NE")
      CALL section_vals_val_get(section, "NETYPE", explicit=explicit)
      ! CALL section_vals_get(section,"NETYPE",explicit=explicit, error=error) !NOTE CAN ONLY DO 1 SECTION - no repeats of list command!
      IF (explicit) THEN
         ALLOCATE (NEtype, stat=stat)
         CPASSERT(stat == 0)
         !CPPostcondition(stat==0,cp_failure_level,routineP,failure)
         NEtype%id_nr = 1
         NEtype%ref_count = 1
         !NEtype%section => section
         !NEtype%label = label
         !region = do_region_global

         simpar%do_NEforce = .TRUE.
         CALL section_vals_val_get(section, "NETYPE", i_val=NEtype%type_of_NE)

         CALL section_vals_val_get(section, "FX", r_val=NEtype%FX)
         CALL section_vals_val_get(section, "FY", r_val=NEtype%FY)
         CALL section_vals_val_get(section, "FZ", r_val=NEtype%FZ)
! Can't really untangle how the particles are defined - it seems odd that particles are tied into force.
! subsys for our cases will be general case, but seriously guys!
         CALL force_env_get(force_env, subsys=subsys)
         CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles)
         atomic_kind_set => atomic_kinds%els
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  shell_present=shell_present) !JREDIT
         IF (shell_present) THEN
            CALL cp_abort(__LOCATION__, &
                          "NETYPE not compatible with shell_present: particles only!")
            !write(*,*) 'Shell not implimented with NE'
            !stop Should be handled by cp_assert!!!!
         ENDIF

!Trying to import thermal_region code to produce particle charges!  Note, not sure how we got derived type attached to particle?!?
!Will assign charges to particles
         particles%els(:)%ccharge = 0
         CALL section_vals_val_get(section, "NELIST", &
                                   i_vals=tmplist)

         IF (NEtype%type_of_NE .EQ. do_ne_colour .AND. MOD(SIZE(tmplist), 2) .NE. 0) THEN
            CALL cp_abort(__LOCATION__, "Colour field has uneven number of partilces!!")
         ENDIF
         !CALL reallocate(t_region%part_index,1,t_region%npart+SIZE(tmplist))
         DO i = 1, SIZE(tmplist)
            ipart = tmplist(i)
            ! TODO fix this to follow current conventions
            !CPPostcondition(((ipart>0).AND.(ipart<=particles%n_els)),cp_failure_level,routineP,error,failure)
            IF (NEtype%type_of_NE .EQ. do_ne_colour) THEN
               particles%els(ipart)%ccharge = (-1)**(i)
            ELSEIF (NEtype%type_of_NE .EQ. do_ne_field) THEN
               particles%els(ipart)%ccharge = 1
            ENDIF
         END DO

         IF (NEtype%type_of_NE .EQ. do_ne_colour .AND. SUM(particles%els(:)%ccharge) .NE. 0) THEN
            CALL cp_abort(__LOCATION__, "Colour charge not balanced!!")
         ENDIF
         CALL NE_info(NETYPE, particles, section)

      ENDIF
   END SUBROUTINE create_NE_type

! TODO: do we need a release_NE_type subroutine to free resources?
! TODO: update the doxygen headers
! *****************************************************************************
!> \brief ...
!> \param nevar ...
!> \param particles ...
!> \param section ...
!> \par History
!>      10.2007 created [tlaino] 06/2017 JREID
!> \author Teodoro Laino, modified James Reid
! **************************************************************************************************
   SUBROUTINE NE_info(nevar, particles, section)
      TYPE(NE_type), POINTER                             :: nevar
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(section_vals_type), POINTER                   :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'NE_info', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iw
      TYPE(cp_logger_type), POINTER                      :: logger

      !failure = .FALSE.
      !NULLIFY(logger)
      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, section, "PRINT%NE_INFO", extension=".log")
      IF (iw > 0) THEN
         WRITE (iw, '( /,A )') ' NON-EQUILIBRIUM| NE Info for MD'
         IF (nevar%type_of_NE .EQ. do_ne_colour) THEN
            WRITE (iw, '( A,T69,A )') ' NON-EQUILIBRIUM| Type of thermostat', 'Colour Field'
         ELSE
            WRITE (iw, '( A,T65,A )') ' NON-EQUILIBRIUM| Type of thermostat', 'Poiseuille Field'
         ENDIF
         WRITE (iw, '( A,A5,A,T69,F12.9 )') &
            ' NON-EQUILIBRIUM| Force in X Direction [', 'au_f', '] ', nevar%FX
         WRITE (iw, '( A,A5,A,T69,F12.9 )') &
            ' NON-EQUILIBRIUM| Force in Y Direction [', 'au_f', '] ', nevar%FY
         WRITE (iw, '( A,A5,A,T69,F12.9 )') &
            ' NON-EQUILIBRIUM| Force in Z Direction [', 'au_f', '] ', nevar%FZ

         WRITE (iw, '(A,/)') ' NON-EQUILIBRIUM| NE "Charges":'

         WRITE (iw, '( A,T72,A,/ )') ' Atom', 'NE charge'
         DO i = 1, SIZE(particles%els(:))
            WRITE (iw, '( I6,T79,I2 )') i, particles%els(i)%ccharge
         ENDDO
         WRITE (iw, '(/,/)')
      END IF
      CALL cp_print_key_finished_output(iw, logger, section, "PRINT%NE_INFO")

   END SUBROUTINE NE_info

END MODULE nonequilibrium_types
