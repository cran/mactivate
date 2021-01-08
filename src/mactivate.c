


#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>








void dmss_dW_b(int *nin, int *din, int *min, double *xMx, double *W, double *MXstar, double *yerrs, double *cc, double *Wout)
{
    int n = nin[0] ;
    int d = din[0] ;
    int m = min[0] ;
    int i ;
    int jm ;
    int jd ;
    int jout ;
    
    int jout_n_key ;
    int jm_n_key ;
    //int j_key ;
    
    double wthis_w1[d] ;
    // long double cumprod ;
    long double d_term ;
    long double f_term ;
    long double BIG_SUM ;
    // W is d x m
    //
    
    for(jout=0; jout<d; jout++) { // big outer loop over columns of X
        
        jout_n_key = jout * n ;
        
        for(jm=0; jm<m; jm++) { // across columns of W
            
            jm_n_key = jm * n ;
            
            for(jd=0; jd<d; jd++) { // down rows of W
                wthis_w1[jd] = W[ jm * d + jd ] ;
            }
            
            BIG_SUM = 0 ;
            for(i=0; i<n; i++) {
                
                d_term = ( xMx[ jout_n_key + i ] - 1 ) ;
                f_term = ( ( xMx[ jout_n_key + i ] * wthis_w1[ jout ] ) + ( 1 - wthis_w1[ jout ] ) ) ;
                
                BIG_SUM = BIG_SUM + MXstar[ jm_n_key + i ] * yerrs[ i ] * d_term / f_term ;
                
            }
            
            Wout[ jm * d + jout ] = 2 * BIG_SUM * cc[ jm ] ;
            
        }
        
        // now MXstar is ready for this particular column of X
        
    }
    
}





void dmss_dW_a(int *nin, int *din, int *min, double *xMx, double *W, double *MXstar, double *yerrs, double *cc, double *Wout)
{
    int n = nin[0] ;
    int d = din[0] ;
    int m = min[0] ;
    int i ;
    int jm ;
    int jd ;
    int jout ;
    //int j_key ;
    
    double wthis_w1[d] ;
    long double cumprod ;
    long double BIG_SUM ;
    // W is d x m
    //
    
    for(jout=0; jout<d; jout++) { // big outter loop over columns of X
        
        for(jm=0; jm<m; jm++) { // across columns of W
            
            for(jd=0; jd<d; jd++) { // down rows of W
                wthis_w1[jd] = W[ jm * d + jd ] ;
            }
            
            
            for(i=0; i<n; i++) {
                
                cumprod = 1 ;
                for(jd=0; jd<d; jd++) {
                    
                    if(jd == jout) {
                        cumprod = cumprod * ( xMx[ jd * n + i ] - 1 ) ;
                    } else {
                        cumprod = cumprod * ( ( xMx[ jd * n + i ] * wthis_w1[ jd ] ) + ( 1 - wthis_w1[ jd ] ) ) ;
                    }
                    
                    
                }
                
                MXstar[ jm * n + i ] = cumprod ;
                
            }
            
        }
        
        // now MXstar is ready for this particular column of X
        
        
        for(jm=0; jm<m; jm++) { // across columns of W
            
            BIG_SUM = 0 ;
            for(i=0; i<n; i++) {
                
                BIG_SUM = BIG_SUM + ( MXstar[ jm * n + i ] * yerrs[ i ] ) ;
                
            }
            
            Wout[ jm * d + jout ] = 2 * BIG_SUM * cc[ jm ] ;
            
        }
        
    }
    
}





void mactivate_a(int *nin, int *din, int *min, double *xMx, double *W, double *MXOUT)
{
    int n = nin[0] ;
    int d = din[0] ;
    int m = min[0] ;
    int i ;
    int jm ;
    int jd ;
    //int j_key ;
    
    double wthis_w1[d] ;
    long double cumprod ;
    
    // W is d x m
    
    for(jm=0; jm<m; jm++) { // across columns of W
        
        for(jd=0; jd<d; jd++) { // down rows of W
            wthis_w1[jd] = W[ jm * d + jd ] ;
        }
        
        
        for(i=0; i<n; i++) {
            
            cumprod = 1 ;
            for(jd=0; jd<d; jd++) {
                cumprod = cumprod * ( ( xMx[ jd * n + i ] * wthis_w1[ jd ] ) + ( 1 - wthis_w1[ jd ] ) ) ;
            }
            
            MXOUT[ jm * n + i ] = cumprod ;
            
        }
        
    }
}







