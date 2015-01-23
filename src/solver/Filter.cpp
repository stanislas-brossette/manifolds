#include <iostream>
#include <pgsolver/solver/Filter.h>
//#include "lex_assert.h"

using namespace Eigen;

namespace pgs
{
  Filter::Filter()
  {
  }

  Filter::Filter(double gamma, Filter::eOption opt)
    : gamma(gamma)
    , option(opt)
  {
  }

  bool Filter::accepts(const VectorXd& p) const
  {
    Entry e = std::make_pair(p,p.norm());
    for (std::list<Entry>::const_iterator it = filter.begin(); it!= filter.end(); ++it)
    {
      if (!accepts(*it,e))
        return false;
    }
    return true;
  }

  bool Filter::dominates(const VectorXd& p) const
  {
    return !accepts(p);
  }

  void Filter::add(const VectorXd& p)
  {
    //lex_assert(accepts(p));
    assert(accepts(p));
    Entry e = std::make_pair(p,p.norm());
    for (std::list<Entry>::iterator it = filter.begin(); it!= filter.end();)
    {
      if (isStronglyDominated(*it, e))
        it = filter.erase(it);
      else
        ++it;
    }
    filter.push_back(e);
  }

  double Filter::getGamma() const
  {
    return gamma;
  }
  void Filter::setGamma(const double newGamma)
  {
    gamma = newGamma;
  }

  Filter::eOption Filter::getOption() const
  {
    return option;
  }
  void Filter::setOption(const Filter::eOption& newOption)
  {
    option = newOption;
  }

  size_t Filter::size() const
  {
    return filter.size();
  }

  const VectorXd& Filter::get(size_t i) const
  {
    std::list<Entry>::const_iterator it = filter.begin();
    advance(it,i);
    return it->first;
  }


  bool Filter::accepts(const Entry& i, const Entry& p) const
  {
    double viol;

    switch(option)
    {
      case EXISTING:  viol = i.second;                    break;
      case TRIAL:     viol = p.second;                    break;
      case MIN:       viol = std::min(i.second,p.second); break;
      case SEPARATE:  return (p.first.array()<(1-gamma)*i.first.array()).any();break;
    }

    return ((i.first - p.first).array()>gamma*viol).any();  // eq. (2.1)
  }

  bool Filter::isStronglyDominated(const Entry& i, const Entry& p) const
  {
    //implements eq. (2.5)
    double a;
    switch(option)
    {
      case EXISTING:  a = gamma * (i.second-p.second);  break;
      case TRIAL:     a = 0;                            break;
      case MIN:       a = gamma * i.second;             break;
      case SEPARATE:  return (i.first.array()>=p.first.array()).all(); break;
    }
    return ((i.first - p.first).array() >= a).all();
  }

  void Filter::print() const
  {
    if (this->size() == 0)
      std::cout << "--empty--" << std::endl;
    for (std::list<Entry>::const_iterator it = filter.begin(); it!= filter.end(); ++it)
    {
      std::cout << it->first.transpose() << std::endl;
    }
  }


  void testFilter01()
  {
    Filter f(1e-5, Filter::EXISTING);
    Vector2d v1(1, 1);
    Vector2d v2(1.1, 1);
    Vector2d v3(1.4, 0.5);
    Vector2d v4(0.8, 1.2);
    Vector2d v5(0.9, 0.9);
    Vector2d v6(0.7, 0.6);

    f.accepts(v1);
    f.add(v1);
    f.accepts(v2);
    f.accepts(v3);
    f.add(v3);
    f.add(v4);
    f.add(v5);
    f.add(v6);
  }
}
