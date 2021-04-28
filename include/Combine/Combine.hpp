/**
 * @file Combine.hpp
 * @author Yifan Zhang (scirocco_gti@yeah.net)
 * @brief
 * @date 2020-11-28 15:46:41
 * @modified: 2021-04-27 15:11:16
 */

#include <bitset>
#include <vector>

std::bitset<64> flip_bit(std::bitset<64> I, uint8_t bit) {
    I[bit] = 1 - I[bit];
    std::bitset<64> ret(I);
    return ret;
}

/**
 * @brief Combine(N, R)
 *
 */
class Combine {
  private:
    const int N = 64;
    int R;
    std::vector<int> pointers;
    int loop(int p);

  public:
    Combine(int);
    ~Combine();
    u_int64_t Next();
    u_int64_t Fetch();
};

Combine::Combine(int r) {
    R = r;
    // init to 0, 1, 2, ..., r-1
    for (int i = 0; i < r; i++) {
        pointers.push_back(i);
    }
}

Combine::~Combine() {}

u_int64_t Combine::Fetch() {
    std::bitset<64> I(0);  // N = 64
    for (int i = 0; i < R; i++) {
        I = flip_bit(I, pointers[i]);
    }
    return I.to_ullong();
}

/**
 * @brief fetch next uint64
 *
 * @return u_int64_t : return 0 if reaches the end
 */
u_int64_t Combine::Next() {
    int loopRet = loop(R - 1);
    if (loopRet == -1) {  // reaches the end
        return 0;
    } else {
        return Fetch();
    }
}

/**
 * @brief recursively move pointers
 *
 * @param p
 * @return int: return 0 if ok, -1 if reached the end
 */
int Combine::loop(int p) {
    pointers[p]++;
    if (pointers[p] > N + p - R) {
        if (p >= 1) {
            if (!loop(p - 1)) {
                for (int i = p; i < R; i++) {
                    pointers[i] = pointers[i - 1] + 1;
                }
                return 0;
            } else {
                return -1;
            }
        } else {  // reached the end
            return -1;
        }
    }
    return 0;
}
