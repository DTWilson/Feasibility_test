// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <Rcpp.h>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <RcppNumerical.h>

using namespace Numer;
using boost::math::negative_binomial_distribution;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_comp_ocs_cpp3(NumericMatrix prm, double n, double mu, int n_e, int n_t)
{
  // Here, calculating power of the main trial as our statistic
  int rows = prm.nrow();
  
  NumericVector results(rows);
  
  for(int i=0; i<rows; i++){
    double total = 0;
    NumericVector sol = prm(i, _);
    double c = sol(0);
    NumericVector v = NumericVector::create(1,2,3);
    NumericVector pr = sol[v];
    // Rcout << pr << std::endl;
    // pr = (p_r, p_a|f, p_f)
    
    // Calculate and store the expected total sample size of 
    // the main trial based on the observed s in the pilot,
    // where s is the number of non-consentors
    double r_max = R::qnbinom(0.999, 2*n, pr(0), 1, 0);
    NumericVector exp_n(r_max+1); 
    NumericVector p_rs(r_max+1); 
    for(int s=r_max; s>=0; s--){
      double r_est = 2*n/(s+2*n);
      double x = 0;
      for(int m=0; m<n_t; m++){
        x = x + R::dbinom(m, n_e, r_est, 0)*m;
      }
      //if(s = 50) Rcout << x << std::endl;
      exp_n(s) = x + n_t*(1 - R::pbinom(n_t-1, n_e, r_est, 1, 0));
      p_rs(s) =  R::dnbinom(s, 2*n, pr(0), 0);
    }
    
    double sd = 1;
    double p_a, p_f, p_r;
    
    for(int a=0; a<=n; a++){
      p_a = R::dbinom(a, n, pr(1), 0); 
      double a_est = a/n;
      
      for(int f=0; f<=2*n; f++){
        p_f = R::dbinom(f, 2*n, pr(2), 0);
        double f_est = f/(2*n);
        
        for(int s=r_max; s>=0; s--){
          
          p_r = p_rs(s); 
          double r_est = (2*n/(2*n+s));
          
          double stat = mu*a_est*sqrt(exp_n(s)*f_est)/sqrt((4 + 2*mu*mu*a_est*(1-a_est)));
          
          //if(r_est*f_est*a_est*a_est > c){
          if(stat > c){
            total = total + p_r*p_a*p_f;
          }
        }
      }
    }
    results(i) = total;
  }
  return results;
}



// [[Rcpp::export]]
NumericVector get_comp_ocs_cpp5(NumericMatrix prm, NumericVector exp_ns, double n, double mu, int n_e, int n_t)
{
  // Here, calculating power of the main trial as our statistic
  int rows = prm.nrow();
  double r_max = exp_ns.size() - 1;
  
  NumericVector results(rows);
  
  for(int i=0; i<rows; i++){
    double total = 0;
    NumericVector sol = prm(i, _);
    double c = sol(0);
    NumericVector v = NumericVector::create(1,2,3);
    NumericVector pr = sol[v];
    // Rcout << pr << std::endl;
    // pr = (p_r, p_a|f, p_f)
    
    NumericVector pr_s(r_max+1);
    for(int s=0; s<=r_max; s++){
      pr_s(s) =  R::pnbinom(s, 2*n, pr(0), 1, 0);
    }
    
    double sd = 1;
    double p_a, p_f, p_r;
    
    for(int a=1; a<=n; a++){
      p_a = R::dbinom(a, n, pr(1), 0);
      double a_est = a/n;

      for(int f=1; f<=2*n; f++){
        p_f = R::dbinom(f, 2*n, pr(2), 0);
        double f_est = f/(2*n);

        double n_cut = c*c*(4 + 2*mu*mu*a_est*(1-a_est))/(mu*mu*a_est*a_est*f_est);

        NumericVector sub = pr_s[exp_ns > n_cut];
        
        //Rcout << f << ", " << sub << std::endl;

        double p_r;
        if(sub.size() != 0){
          p_r = sub(sub.size()-1);
          //Rcout << p_r << std::endl;
        } else {
          p_r = 0;
        }

        total = total + p_r*p_a*p_f;
      }
    }
    results(i) = total;
  }
return results;
}


// [[Rcpp::export]]
NumericVector get_comp_ocs_cpp_precomp(double c, NumericMatrix pr_m, NumericMatrix pa_m, NumericMatrix pf_m, 
                                       NumericVector exp_ns,
                                       double n, double mu, int n_e, int n_t)
{
  // Here, calculating power of the main trial as our statistic
  int rows = pr_m.nrow();
  
  NumericVector results(rows);
  
  int r_max = pr_m.ncol()-1;
  
  for(int i=0; i<rows; i++){
    
    double total = 0;
    double sd = 1;
    double p_a, p_f, p_r;
    
    for(int a=0; a<=n; a++){
      p_a = pa_m(i,a+1);
      double a_est = a/n;
      
      for(int f=0; f<=2*n; f++){
        p_f = pf_m(i,f+1);
        double f_est = f/(2*n);
        
        for(int s=r_max; s>0; s--){
          p_r = pr_m(i,s); 
          double r_est = (2*n/(2*n+s));
          Rcout << s << std::endl;
          double stat = mu*a_est*sqrt(exp_ns(s-1)*f_est)/sqrt((4 + 2*mu*mu*a_est*(1-a_est)));
          
          if(stat > c){
            total = total + p_r*p_a*p_f;
          }
        }
      }
    }
    results(i) = total;
  }
  return results;
}

// [[Rcpp::export]]
NumericMatrix vec_f(NumericMatrix sols, double n, double x0, double x1, double mu, int n_e, int n_t)
{
  int rows = sols.nrow();
  NumericMatrix results(2, rows);
  
  // Get the type I errors
  NumericMatrix tI_mat(rows,4);
  tI_mat(_,0) = sols(_,0);
  tI_mat(_,1) = sols(_,1);
  tI_mat(_,2) = sols(_,2);
  tI_mat(_,3) = sols(_,3);
  
  results(0,_) = get_comp_ocs_cpp3(tI_mat, n, mu, n_e, n_t);
  
  // Get the type II errors
  NumericMatrix tII_mat(rows,4);
  tII_mat(_,0) = sols(_,0);
  tII_mat(_,1) = sols(_,4);
  tII_mat(_,2) = sols(_,5);
  tII_mat(_,3) = sols(_,6);

  
  results(1,_) = 1-get_comp_ocs_cpp3(tII_mat, n, mu, n_e, n_t);
  
  return results;
}

// [[Rcpp::export]]
NumericMatrix vec_g(NumericMatrix sols, double n, double x0, double x1, double mu, int n_e, int n_t)
{
  int rows = sols.nrow();
  NumericMatrix results(2, rows);
  
  for(int i=0; i<rows; i++){
    NumericVector sol = sols(i, _);
    double c = sol(0);
    NumericVector vI = NumericVector::create(1,2,3);
    NumericVector prI = sol[vI];
    // pr = (p_r, p_a|f, p_f)
    NumericVector vII = NumericVector::create(4,5,6);
    NumericVector prII = sol[vII];
    
    double exp_nI = n_t*(1-R::pbinom(n_t-1, n_e, prI(0),1, 0));
    for(int j=0; j<n_t; j++){
      exp_nI = exp_nI + R::dbinom(j, n_e, prI(0), 0)*j;
    }
    
    results(0,i) = x0 - prI(1)*mu*sqrt(prI(2)*exp_nI)/sqrt(2 + mu*mu * prI(1)*(1-prI(1)));

    
    double exp_nII = n_t*(1-R::pbinom(n_t-1, n_e, prII(0), 1, 0));
    for(int j=0; j<n_t; j++){
      exp_nII = exp_nII + R::dbinom(j, n_e, prII(0), 0)*j;
    }
    
    results(1,i) = prII(1)*mu*sqrt(prII(2)*exp_nII)/sqrt(2 + mu*mu * prII(1)*(1-prII(1))) - x1;
  }
  
  return results;
}

// [[Rcpp::export]]
NumericVector get_comp_ocs_cpp4(double c, double n, NumericMatrix prm, double mu, int n_e, int n_t)
{
  // Here, calculating power of the main trial as our statistic
  // Allowing for crrelation between adherence and follow-up outcomes
  int rows = prm.nrow();
  double total = 0;
  
  NumericVector n_r(n_t);
  for(int i=0; i<n_t; i++){
    n_r(i) = i;
  }
  
  NumericVector results(rows);
  
  for(int i=0; i<rows; i++){
    NumericVector pr = prm(i, _);
    // pr = (p_r, p_a|f, p_f, OR)
    NumericVector probs(4);
    // probs = (p_af, p_!af, p_a!f, p_!a!f)
    probs(0) = pr(1)*pr(2);
    probs(1) = pr(2)-probs(0);
    probs(3) = probs(1)*pr(3)*(1-probs(0)-probs(1))/(probs(0) + probs(1)*pr(3));
    //probs(3) = (pr(3)*probs(1)*(1-probs(0)-probs(1))/probs(0))/(1+pr(3)*probs(1)/probs(0));
    probs(2) = 1 - probs(0) - probs(1) - probs(3);
    
    double r_max = R::qnbinom(0.999, 2*n, pr(0), 1, 0);
    NumericVector exp_n(r_max+1); 
    NumericVector p_rs(r_max+1); 
    for(int s=r_max; s>=0; s--){
      double r_est = 2*n/(s+2*n);
      double x = 0;
      for(int m=0; m<n_t; m++){
        x = x + R::dbinom(n_r(m), n_e, r_est, 0)*n_r(m);
      }
      exp_n(s) = x + n_t*(1 - R::pbinom(n_t-1, n_e, r_est, 1, 0));
      p_rs(s) =  R::dnbinom(s, 2*n, pr(0), 0);
    }
    
    double sd = 1;
    // go through the possible adherence x follow-up outcomes in the intervention arm
    for(int af=0; af<=n; af++){
      for(int naf=0; naf<=(n-af); naf++){
        for(int anf=0; anf<=(n-af-naf); anf++){
          int nanf = n-af-naf-anf;
          NumericVector afs(4);
          afs(0) = af; afs(1) = naf; afs(2) = anf; afs(3) = nanf;
          double p_int1 = boost::math::factorial<double>(afs(0))*boost::math::factorial<double>(afs(1))*boost::math::factorial<double>(afs(2))*boost::math::factorial<double>(afs(3));
          double p_int2 = pow(probs(0), afs(0))*pow(probs(1), afs(1))*pow(probs(2), afs(2))*pow(probs(3), afs(3));
          double p_int= boost::math::factorial<double>(n)/p_int1*p_int2;
          
          //Rcout << afs << " - " << p_int << std::endl;
          
          double a_est = (af+anf)/n;
          
          // go through the follow-up outcomes in the control arm
          for(int f_c=0; f_c<=n; f_c++){
            double p_fc = R::dbinom(f_c, n, pr(2), 0);
            double f = af+naf+f_c;
            double f_est = f/(2*n);
            
            // go through the recruitment outcomes
            for(int s=r_max; s>=0; s--){
              
              double p_r = p_rs(s); 
              double r_est = (2*n/(2*n+s));
              
              double stat = mu*a_est*sqrt(exp_n(s)*f_est)/sqrt((2 + mu*mu*a_est*(1-a_est)));

              if(stat > c){
                total = total + p_r*p_int*p_fc;
              }
            }
          }
        }
      }
    }
    results(i) = total;
  }
  return results;
}



// [[Rcpp::export]]
NumericVector get_ocs_var_cpp(double c, double n, NumericMatrix prm, double mu, int n_e, int n_t)
{
  // Here, calculating power of the main trial as our statistic
  int rows = prm.nrow();
  double total = 0;
  
  // fill a vector with possible numbers recruited
  NumericVector n_r(n_t);
  for(int i=0; i<n_t; i++){
    n_r(i) = i;
  }
  
  NumericVector results(rows);
  
  for(int i=0; i<rows; i++){
    NumericVector pr = prm(i, _);
    // pr = (p_r, p_a|f, p_f, sd)
    
    // find the exp_n corresponding to each possible value of 
    // the number of eligibles not recruited, so we can look
    // up later
    double r_max = R::qnbinom(0.999, 2*n, pr(0), 1, 0);
    NumericVector exp_n(r_max+1); 
    NumericVector p_rs(r_max+1); 
    for(int s=r_max; s>=0; s--){
      double r_est = 2*n/(s+2*n);
      double x = 0;
      for(int m=0; m<n_t; m++){
        x = x + R::dbinom(n_r(m), n_e, r_est, 0)*n_r(m);
      }
      exp_n(s) = x + n_t*(1 - R::pbinom(n_t-1, n_e, r_est, 1, 0));
      p_rs(s) =  R::dnbinom(s, 2*n, pr(0), 0);
    }
    
    double p_a, p_f, p_r, p_var;
    double a_est, f_est, r_est, var_crit;
    
    for(int a=0; a<=n; a++){
      p_a = R::dbinom(a, n, pr(1), 0); 
      a_est = a/n;
      
      for(int f=1; f<=2*n; f++){
        p_f = R::dbinom(f, 2*n, pr(2), 0);
        f_est = f/(2*n);
        
        for(int s=r_max; s>=0; s--){

          p_r = p_rs(s);
          r_est = (2*n/(2*n+s));

          var_crit = (mu*mu*a_est*a_est*f_est*exp_n(s) - 2*mu*mu*a_est*(1-a_est)*c*c)/(4*c*c);
          p_var = R::pchisq(var_crit*(f-1)/pow(pr(3), 2.0), f-1, 1, 0);

          //Rcout << r_est << ", " << a_est << ", " << f_est << ", " << var_crit << std::endl;
          //Rcout << p_r << ", " << p_a << ", " << p_f << ", " << p_var << std::endl;

          total = total + p_var*p_r*p_a*p_f;
          }
        }
      }
    results(i) = total;
    }
  return results;
}


// [[Rcpp::export]]
double MC_cpp(double MC, double c, NumericVector pr, double n, double mu, int n_e, int n_t)
{
  // pr = (p_r, p_a|f, p_f, sd)
  
  NumericVector n_r(n_t);
  for(int i=0; i<n_t; i++){
    n_r(i) = i;
  }
  
  double r_max = R::qnbinom(0.999, 2*n, pr(0), 1, 0);
  NumericVector exp_n(r_max+1); 
  for(int s=r_max; s>=0; s--){
    double r_est = 2*n/(s+2*n);
    
    double x = 0;
    for(int m=0; m<n_t; m++){
      x = x + R::dbinom(n_r(m), n_e, r_est, 0)*n_r(m);
    }
    
    exp_n(s) = x + n_t*(1 - R::pbinom(n_t-1, n_e, r_est, 1, 0));
  }
 
  double sd = 1;
  
  NumericVector a = Rcpp::rbinom(MC, n, pr(1));
  NumericVector f = Rcpp::rbinom(MC, 2*n, pr(2));
  NumericVector s = clamp(0, Rcpp::rnbinom(MC, 2*n, pr(0)), r_max);
  
  double a_est, f_est;
  
  double total = 0;
  
  for(int i=0; i < MC; i++){
    a_est = a(i)/n;
    f_est = f(i)/(2*n);
    
    double stat = mu*a_est*sqrt(exp_n(s(i))*f_est)/sqrt((4 + 2*mu*mu*a_est*(1-a_est)));
    
    //Rcout << stat << std::endl;
    
    if(stat > c){
      total = total + 1;
    }
  }

  return total/MC;
}




//to_sum <- NULL
//for(i in 2:length(mus)){
//  to_sum <- c(to_sum, 0.5*(mus_probs[i-1] + mus_probs[i])*(mus[i]-mus[i-1]))
//}


// see https://cran.r-project.org/web/packages/RcppNumerical/vignettes/introduction.html

