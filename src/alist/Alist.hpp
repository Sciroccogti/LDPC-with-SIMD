/*
 * @Author: Yifan Zhang
 * @Date: 2020-10-29 15:30:45
 * @Last Modified by: Yifan Zhang
 * @Last Modified time: 2020-10-29 15:43:00
 */

#include "alist_matrix.h"
#include "nbalist_matrix.h"

template <class T>
class Alist {
  private:
    T data;

  public:
    Alist(T d);
    ~Alist();
};
