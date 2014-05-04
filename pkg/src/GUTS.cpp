/** GUTS Class Implementation.
 *
 * @file		GUTS.cpp
 * @author  soeren.vogel@posteo.ch
 * @author	carlo.albert@eawag.ch
 * @date		2012-05-11
 * @license	GPL-2
 */

#ifndef guts_h
#include "GUTS.h"
#endif

using namespace std;

/*
 * Constructor.
 */
GUTS::GUTS() : 
  mTitle("Experiment"),
  mC(), mCt(), my(), myt(),
  mpar(), mM(10000), mdist("lognormal"), mN(10000),
  mz(), mD(), mS(), mLL(GNAN_DOUBLE),
  mErrors(GER_COUNT, true), mErrorMessages(GER_COUNT, "")
{
  mErrorMessages.at(GER_C)    = "C not set up";
  mErrorMessages.at(GER_CT)   = "Ct not set up";
  mErrorMessages.at(GER_Y)    = "y not set up";
  mErrorMessages.at(GER_YT)   = "yt not set up";
  mErrorMessages.at(GER_PAR)  = "par not set up";

  // M has a default value of 10000
  mErrors.at(GER_M)           = false;
  mErrorMessages.at(GER_M)    = "";
  
  // dist has a valid default:
  mErrors.at(GER_DIST)        = false;
  mErrorMessages.at(GER_DIST) = "";
  
  // N has a default value of 10000
  mErrors.at(GER_N)           = false;
  mErrorMessages.at(GER_N)    = "";

  mErrorMessages.at(GER_Z)    = "z not available";
  
  // There's ToDos on D, use no error so far
  mErrors.at(GER_D)           = false;
  mErrorMessages.at(GER_D)    = "";
  
  mErrorMessages.at(GER_S)    = "Survival probabilities not calculated";
  mErrorMessages.at(GER_LL)   = "Loglikelihood not calculated";
} // end GUTS::GUTS()

/*
 * Destructor.
 */
GUTS::~GUTS() {
  mC.clear();
  mCt.clear();
  my.clear();
  myt.clear();
  mpar.clear();
  mz.clear();
  mD.clear();
  mS.clear();
  mErrors.clear();
  mErrorMessages.clear();
} // end GUTS::~GUTS()

/*
 * mC and mCt
 * Error read: none
 * Error write: GER_C, GER_CT
 */
void GUTS::setConcentrations(const vector<double>& C, const vector<double>& Ct) {
  // Reset errors on this method (false means no error).
  int i[] = {GER_C, GER_CT};
  for (unsigned int j = 0; j < 2; ++j) {
    mErrors.at(i[j])        = false;
    mErrorMessages.at(i[j]) = "";
  }

  // New errors on C
  if (C.size() < 2) {
    mErrors.at(GER_C)         = true;
    mErrorMessages.at(GER_C)  = "vector C must be longer than 1";
  }
  if (C.size() != Ct.size()) {
    mErrors.at(GER_C)         = true;
    mErrorMessages.at(GER_C)  = "vectors C and Ct must have the same length";
  }

  // New errors on Ct
  if ( Ct.size() < 2 )
  {
    mErrors.at( GER_CT )        = true;
    mErrorMessages.at( GER_CT ) = "Vector Ct must be longer than 1.";
  }
  if ( Ct.size() != C.size() )
  {
    mErrors.at( GER_CT )        = true;
    mErrorMessages.at( GER_CT ) = "Vectors Ct and C must have the same length.";
  }
  if ( Ct.at(0) != 0.0f )
  {
    mErrors.at( GER_CT )        = true;
    mErrorMessages.at( GER_CT ) = "Vector Ct must start at 0.0.";
  }
  for ( unsigned int i = 1; i < Ct.size(); ++i )
  {
    // Prevent division by 0.
    if ( (Ct.at(i) - Ct.at(i-1)) <= 0.0f )
    {
      mErrors.at( GER_CT )		= true;
      mErrorMessages.at( GER_CT ) = "Vector Ct must contain unique values in ascending order.";
    }
  }
  
  // Assign C if no errors, no errors means false!
  if ( ! mErrors.at( GER_C ) )
    mC = C;
  else
    mC.resize( 0 );
  
  // Assign Ct if no errors, no errors means false!
  if ( ! mErrors.at( GER_CT ) )
    mCt = Ct;
  else
    mCt.resize( 0 );
} // end GUTS::setConcentrations()

/*
 * my and myt.
 * error read: none
 * error write: GER_Y, GER_YT
 */
void GUTS::setSurvivors( const vector<int>& y, const vector<double>& yt )
{
  // Reset errors on this method.
  int i[] = { GER_Y, GER_YT };
  for ( unsigned int j = 0; j < 2; ++j )
  {
    mErrors.at( i[j] )          = false; // false means no error
    mErrorMessages.at( i[j] )   = "";
  }
  
  // New errors on y
  if ( y.size() < 2 )
  {
    mErrors.at( GER_Y )			= true;
    mErrorMessages.at( GER_Y )  = "Vector y must be longer than 1.";
  }
  if ( y.size() != yt.size() )
  {
    mErrors.at( GER_Y )         = true;
    mErrorMessages.at( GER_Y )  = "Vectors y and yt must have the same length.";
  }
  
  // New errors on yt
  if ( yt.size() < 2 )
  {
    mErrors.at( GER_YT )        = true;
    mErrorMessages.at( GER_YT ) = "Vector yt must be longer than 1.";
  }
  if ( yt.size() != y.size() )
  {
    mErrors.at( GER_YT )        = true;
    mErrorMessages.at( GER_YT ) = "Vectors yt and y must have the same length.";
  }
  if ( yt.at(0) != 0.0f )
  {
    mErrors.at( GER_YT )        = true;
    mErrorMessages.at( GER_YT ) = "Vector yt must start at 0.0.";
  }
  for ( unsigned int i = 1; i < yt.size(); ++i )
  {
    // Prevent division by 0
    if ( (yt.at(i) - yt.at(i-1)) <= 0.0f )
    {
      mErrors.at(GER_YT)			= true;
      mErrorMessages.at(GER_YT)  	= "Vector yt must contain unique values in ascending order.";
    }
  }
  
  // Assign y if no errors
  if ( ! mErrors.at( GER_Y ) )
    my = y;
  else
    my.resize( 0 );
  
  // Assign yt if no errors
  if ( ! mErrors.at( GER_YT ) )
    myt = yt;
  else
    myt.resize( 0 );
} // end GUTS::setSurvivors()

/*
 * mpar.
 * error read: none
 * error write: GER_PAR
 * Sample parameters will be checked in calcSample()!
 */
void GUTS::setParameters( const vector<double>& par )
{
  // Reset error on this method.
  mErrors.at( GER_PAR )          = false; // false means no error
	mErrorMessages.at( GER_PAR )   = "";
  
  // New errors on par
  if ( par.size() < 3 )
  {
    mErrors.at( GER_PAR )			= true;
    mErrorMessages.at( GER_PAR )	= "Vector par must be longer than 2.";
  }
  if ( *min_element( par.begin(), par.begin()+3 ) < 0.0 )
  {
    mErrors.at( GER_PAR )			= true;
    mErrorMessages.at( GER_PAR )   	= "par 1, 2 and 3 must contain non-negative values.";
  }
  
  // Assign par if no errors
  if ( ! mErrors.at( GER_PAR ) )
    mpar = par;
  else
    mpar.resize( 0 );
} // end GUTS::setParameters()

/*
 * mM.
 * error read: none
 * error write: GER_M
 */
void GUTS::setTimeGridPoints( const int& M )
{
  // Reset errors on this method.
  mErrors.at( GER_M )			= false; // false means no error
  mErrorMessages.at( GER_M ) 	= "";
  
  // New errors on M
  if ( M < 2 )
  {
    mErrors.at( GER_M )			= true;
    mErrorMessages.at( GER_M ) 	= "M must be greater than 1.";
  }
  
  // Assign M if no errors
  if ( ! mErrors.at( GER_M ) )
    mM = M;
  else
    mM = GNAN_INT;
} // end GUTS::setTimeGridPoints()

/*
 * dist.
 * error read: none
 * error write: GER_DIST
 */
void GUTS::setDistribution( const string dist )
{
  // Reset error on dist
  mErrors.at( GER_DIST )			= false;
  mErrorMessages.at( GER_DIST )   = "";
  
  // New error on dist
  if ( dist != "lognormal" && dist != "empirical" )
  {
    mErrors.at( GER_DIST )			= true;
    mErrorMessages.at( GER_DIST )	= "Unknown distribution name.";
  }
  
  // Assign dist if no errors
  if ( ! mErrors.at( GER_DIST ) )
    mdist = dist;
  else
    mdist = GNAN_STRING;
} // end GUTS::setDistribution()

/*
 * mN.
 * error read: none
 * error write: GER_N, GER_Z
 */
void GUTS::setSampleLength( const int& N )
{
  // Reset errors on this method.
  mErrors.at( GER_N )			= false;
  mErrorMessages.at( GER_N )	= "";
  
  // New errors on N
  if ( N < 3 )
  {
    mErrors.at( GER_N )			= true;
    mErrorMessages.at( GER_N ) 	= "N must be greater than 2";
  }
  
  // Assign n and size of z if no errors
  if ( ! mErrors.at( GER_N ) )
  {
    mN = N;
    mz.resize( mN );
  }
  else
  {
    mErrors.at( GER_Z ) = true;
    mErrorMessages.at( GER_Z ) = "failed to set sample length N";
    mN = GNAN_INT;
  }
} // end GUTS::setSampleLength()

/*
 * mz.
 * error read: none
 * error write: GER_Z
 */
void GUTS::setSample( const vector<double>& z )
{
  // Reset error on z
  mErrors.at( GER_Z )			= false;
  mErrorMessages.at( GER_Z )  = "";
  
  // New error on z
  if ( z.size() < 3 )
  {
    mErrors.at( GER_Z )			= true;
    mErrorMessages.at( GER_Z ) 	= "z must be longer than 2";
  }
  else if ( *min_element( z.begin(), z.end() ) < 0.0f )
  {
    mErrors.at( GER_Z )     	= true;
    mErrorMessages.at( GER_Z )	= "z must contain non-negative values";
  }
  
  // Assign z, dist and N if no errors
  if ( ! mErrors.at( GER_Z ) )
  {
    mz = z;
    setDistribution( "empirical" );
    setSampleLength( mz.size() );
  }
} // end GUTS::setSample()
