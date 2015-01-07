#ifndef _FILTER_H_
#define _FILTER_H_

#include<Eigen/Core>
#include<list>

class Filter
{
public:
  /*Option refer to the different acceptability conditions (2.1) with (2.2) - (2.4)
   *in the paper "A multidimensional filter algorithm for nonlinear equations and 
   *nonlinear least-squares" by N.I.M Gould, S. Leyffer and Ph. L. Toint
   */
  enum eOption
  {
    EXISTING, //2.2
    TRIAL,    //2.3
    MIN,      //2.4
    SEPARATE  //implements the initial u_i<(1-gamma)*v_i of Flechter
  };

  typedef std::pair<Eigen::VectorXd, double> Entry;

  Filter(double gamma, eOption opt=EXISTING);

  bool accepts(const Eigen::VectorXd& p) const;
  bool dominates(const Eigen::VectorXd& p) const;
  void add(const Eigen::VectorXd& p);

  double getGamma() const;
  eOption getOption() const;

  size_t size() const;
  const Eigen::VectorXd& get(size_t i) const;


private:
  bool accepts(const Entry& i, const Entry& p) const;             //i accepts p 
  bool isStronglyDominated(const Entry& i, const Entry& p) const; //i is strongly dominated by p

private:
  eOption option;
  double  gamma;
  std::list<Entry> filter;
};

void testFilter01();

#endif //_FILTER_H_