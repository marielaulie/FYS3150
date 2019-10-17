TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas

  SUBROUTINE gauleg(x1,x2,x,w,n)
    IMPLICIT NONE
    INTEGER :: i, j, m, n
    REAL(DP) :: eps, x1, x2, x, w
    DIMENSION :: x(n), w(n)
    PARAMETER (eps=3.D-14)
    REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.
          p2=0.
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.*xl/((1.-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO

  END SUBROUTINE gauleg

  !     Function to integrate a function func over the
  !     interval [a,b] with input a, b, and the number of steps
  !     n.  It returns the sum as the variable trapez_sum
  !     The trapezoidal rule is used


  SUBROUTINE trapezoidal_rule(a,b,trapez_sum,n,func)
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: a,b
    REAL(DP), INTENT(INOUT) :: trapez_sum
    REAL(DP) fa, fb, x, step
    INTEGER :: j
    INTERFACE
       DOUBLE PRECISION FUNCTION  func(x)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: x

       END FUNCTION func
    END INTERFACE


    step=(b-a)/FLOAT(n)
    fa=func(a)/2. ; fb=func(b)/2. ; trapez_sum=0.
    DO j=1,n-1
       x=j*step+a
       trapez_sum=trapez_sum+func(x)
    ENDDO
    trapez_sum=(trapez_sum+fb+fa)*step

  END  SUBROUTINE trapezoidal_rule


  !     Function to integrate a function func over the
  !     interval [a,b] with input a, b, and the number of steps
  !     n.  It returns the sum as the variable simpson_sum
  !     Simpson's method is used


  SUBROUTINE simpson(a,b,simpson_sum,n,func)
    REAL(DP), INTENT(IN) :: a,b
    REAL(DP), INTENT(INOUT) :: simpson_sum
    REAL(DP) fa, fb, x, step, fac
    INTEGER, INTENT(IN) :: n
    INTEGER :: j

    INTERFACE
       DOUBLE PRECISION FUNCTION  func(x)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: x

       END FUNCTION func
    END INTERFACE

    step=(b-a)/FLOAT(n)
    fa=func(a) ; fb=func(b) ; simpson_sum=fa ; fac=2.
    DO j=1,n-1
       IF ( fac == 2.) THEN
          fac = 4.
       ELSE
          fac = 2.
       ENDIF
       x=j*step+a
       simpson_sum=simpson_sum+func(x)*fac
    ENDDO
    simpson_sum=(simpson_sum+fb)*step/3.

  END SUBROUTINE simpson
