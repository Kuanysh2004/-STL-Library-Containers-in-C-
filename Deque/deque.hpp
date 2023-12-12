#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

template <typename T, typename Allocator = std::allocator<T>>
class Deque {
 public:
  using value_type = std::remove_cv_t<T>;
  using allocator_type = Allocator;
  using alloc_traits = std::allocator_traits<Allocator>;
  Deque() {
    size_ = 0;
    end_ind_1_ = 0;
    end_ind_2_ = 0;
    bedin_ind_1_ = 0;
    begin_ind_2_ = 0;
  }

  Deque(const Allocator& alloc) {
    size_ = 0;
    end_ind_1_ = 0;
    end_ind_2_ = 0;
    bedin_ind_1_ = 0;
    begin_ind_2_ = 0;
    alloc_ = alloc;
  }

  Deque(size_t count, const Allocator& alloc = Allocator()) {
    size_ = count;
    bedin_ind_1_ = 0;
    begin_ind_2_ = 0;
    alloc_ = alloc;
    if (count == 0) {
      end_ind_1_ = 0;
      end_ind_2_ = 0;
    } else {
      end_ind_1_ = (count - 1) / kColumn;
      end_ind_2_ = (count - 1) % kColumn;
    }
    size_t iii = 0;
    size_t jjj = 0;
    size_t constructed = 0;
    try {
      deq_.resize((count + kColumn - 1) / kColumn);
      for (; iii < deq_.size() - 1; ++iii) {
        deq_[iii] = alloc_traits::allocate(alloc_, kColumn);
        for (size_t j; j < kColumn; ++j) {
          alloc_traits::construct(alloc_, deq_[iii] + j);
          ++constructed;
        }
      }
      deq_[iii] = alloc_traits::allocate(alloc_, kColumn);
      for (; jjj < end_ind_2_ + 1; ++jjj) {
        alloc_traits::construct(alloc_, deq_[iii] + jjj);
        ++constructed;
      }
    } catch (...) {
      while (constructed != 0) {
        pop_back();
        --constructed;
      }
      for (size_t i = 0; i < deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kColumn);
      }
      throw;
    }
  }

  Deque(size_t count, const T& value, const Allocator& alloc = Allocator()) {
    alloc_ = alloc;
    size_ = count;
    bedin_ind_1_ = 0;
    begin_ind_2_ = 0;
    if (count == 0) {
      end_ind_1_ = 0;
      end_ind_2_ = 0;
    } else {
      end_ind_1_ = (count - 1) / kColumn;
      end_ind_2_ = (count - 1) % kColumn;
    }
    size_t iii = 0;
    size_t jjj = 0;
    try {
      deq_.resize((count + kColumn - 1) / kColumn);
      for (; iii < deq_.size() - 1; ++iii) {
        deq_[iii] = alloc_traits::allocate(alloc_, kColumn);
        for (size_t j = 0; j < kColumn; ++j) {
          alloc_traits::construct(alloc_, deq_[iii] + j, value);
        }
      }
      deq_[iii] = alloc_traits::allocate(alloc_, kColumn);
      for (; jjj < end_ind_2_ + 1; ++jjj) {
        alloc_traits::construct(alloc_, deq_[iii] + jjj, value);
      }
    } catch (...) {
      while (size_ != 0) {
        pop_back();
      }
      for (size_t i = 0; i < deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kColumn);
      }
      throw;
    }
  }

  ~Deque() {
    if (!empty()) {
      while (size_ != 0) {
        pop_back();
      }
    }
    for (size_t i = 0; i < deq_.size(); ++i) {
      alloc_traits::deallocate(alloc_, deq_[i], kColumn);
    }
  }

  Deque(const Deque& other) {
    alloc_ = alloc_traits::select_on_container_copy_construction(other.alloc_);
    auto iter = other.begin();
    try {
      size_t iii = other.bedin_ind_1_;
      size_t jjj = other.begin_ind_2_;
      deq_.resize(other.deq_.size());
      for (size_t i = 0; i < other.deq_.size(); ++i) {
        deq_[i] = alloc_traits::allocate(alloc_, kColumn);
      }
      for (; iter < other.end(); ++iter) {
        alloc_traits::construct(alloc_, deq_[iii / kColumn] + (jjj % kColumn),
                                *iter);
        ++iii;
        ++jjj;
      }
    } catch (...) {
      size_t iii = other.bedin_ind_1_;
      size_t jjj = other.begin_ind_2_;
      for (auto it = other.begin(); it < iter; ++it) {
        alloc_traits::destroy(alloc_, deq_[iii / kColumn] + (jjj % kColumn));
        ++iii;
        ++jjj;
      }
      for (size_t i = other.bedin_ind_1_;
           i < other.deq_.size() + other.bedin_ind_1_; ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kColumn);
      }
      throw;
    }
    size_ = other.size_;
    bedin_ind_1_ = other.bedin_ind_1_;
    begin_ind_2_ = other.begin_ind_2_;
    end_ind_1_ = other.end_ind_1_;
    end_ind_2_ = other.end_ind_2_;
  }

  Deque& operator=(const Deque<T, Allocator>& other) {
    std::vector<T*> ddd;
    auto iter = other.begin();
    size_t iii = other.bedin_ind_1_;
    size_t jjj = other.begin_ind_2_;
    auto tempalloc = alloc_;
    if (std::allocator_traits<decltype(
            other.alloc_)>::propagate_on_container_copy_assignment::value) {
      alloc_ = other.alloc_;
    }
    try {
      ddd.resize(other.deq_.size());
      for (size_t i = 0; i < other.deq_.size(); ++i) {
        ddd[i] = alloc_traits::allocate(alloc_, kColumn);
      }
      for (; iter < other.end(); ++iter) {
        alloc_traits::construct(alloc_, ddd[iii / kColumn] + (jjj % kColumn),
                                *iter);
        ++iii;
        ++jjj;
      }
    } catch (...) {
      size_t ii2 = other.bedin_ind_1_;
      size_t jj2 = other.begin_ind_2_;
      for (auto it = other.begin(); it < iter; ++it) {
        alloc_traits::destroy(alloc_, ddd[ii2 / kColumn] + (jj2 % kColumn));
        ++ii2;
        ++jj2;
      }
      for (size_t i = other.bedin_ind_1_; i < other.deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, ddd[i], kColumn);
      }
      throw;
    }
    auto tmpalloc = alloc_;
    alloc_ = tempalloc;
    while (size_ != 0) {
      pop_back();
    }
    for (size_t i = 0; i < deq_.size(); ++i) {
      alloc_traits::deallocate(alloc_, deq_[i], kColumn);
    }
    deq_ = ddd;
    alloc_ = tmpalloc;
    size_ = other.size_;
    bedin_ind_1_ = other.bedin_ind_1_;
    begin_ind_2_ = other.begin_ind_2_;
    end_ind_1_ = other.end_ind_1_;
    end_ind_2_ = other.end_ind_2_;
    return *this;
  }

  Deque(Deque&& other) {
    deq_ = std::move(other.deq_);
    alloc_ = std::move(other.alloc_);
    size_ = std::move(other.size_);
    bedin_ind_1_ = std::move(other.bedin_ind_1_);
    begin_ind_2_ = std::move(other.begin_ind_2_);
    end_ind_1_ = std::move(other.end_ind_1_);
    end_ind_2_ = std::move(other.end_ind_2_);
    other.deq_.resize(0);
    other.size_ = 0;
    other.bedin_ind_1_ = 0;
    other.begin_ind_2_ = 0;
    other.end_ind_1_ = 0;
    other.end_ind_2_ = 0;
  }

  Deque& operator=(Deque&& other) {
    auto tmp = std::move(other);
    std::swap(deq_, tmp.deq_);
    std::swap(alloc_, tmp.alloc_);
    std::swap(size_, tmp.size_);
    std::swap(bedin_ind_1_, tmp.nachalo_ind_1_);
    std::swap(begin_ind_2_, tmp.nachalo_ind_2_);
    std::swap(end_ind_1_, tmp.konec_ind_1_);
    std::swap(end_ind_2_, tmp.konec_ind_2_);
    return *this;
  }

  Deque(std::initializer_list<T> init, const Allocator& alloc = Allocator()) {
    alloc_ = alloc;
    auto itt = init.begin();
    deq_.resize((init.size() + kColumn - 1) / kColumn);
    try {
      for (size_t i = 0; i < deq_.size(); ++i) {
        deq_[i] = alloc_traits::allocate(alloc_, kColumn);
      }
      size_t iii = 0;
      size_t jjj = 0;
      for (; itt != init.end(); ++itt) {
        alloc_traits::construct(alloc_, deq_[iii / kColumn] + (jjj % kColumn),
                                *itt);
        ++iii;
        ++jjj;
      }
    } catch (...) {
      size_t iii = 0;
      size_t jjj = 0;
      for (auto it = init.begin(); it != itt; ++it) {
        alloc_traits::destroy(alloc_, deq_[iii / kColumn] + (jjj % kColumn));
        ++iii;
        ++jjj;
      }
      for (size_t i = 0; i < deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kColumn);
      }
      throw;
    }
    size_ = init.size();
    bedin_ind_1_ = 0;
    begin_ind_2_ = 0;
    if (init.size() == 0) {
      end_ind_1_ = 0;
      end_ind_2_ = 0;
    } else {
      end_ind_1_ = (init.size() - 1) / kColumn;
      end_ind_2_ = (init.size() - 1) % kColumn;
    }
  }

  size_t size() const { return size_; }

  bool empty() const { return (size_ == 0); }

  T& operator[](size_t iii) {
    return deq_[(bedin_ind_1_ * kColumn + iii + begin_ind_2_) / kColumn]
               [(begin_ind_2_ + iii) % kColumn];
  }

  T& at(size_t iii) {
    if (iii >= size_) {
      throw std::out_of_range("");
    }
    return deq_[(bedin_ind_1_ * kColumn + iii + begin_ind_2_) / kColumn]
               [(begin_ind_2_ + iii) % kColumn];
  }

  const T& operator[](size_t iii) const {
    return deq_[(bedin_ind_1_ * kColumn + iii + begin_ind_2_) / kColumn]
               [(begin_ind_2_ + iii) % kColumn];
  }

  const T& at(size_t iii) const {
    if (iii >= size_) {
      throw std::out_of_range("");
    }
    return deq_[(bedin_ind_1_ * kColumn + iii + begin_ind_2_) / kColumn]
               [(begin_ind_2_ + iii) % kColumn];
  }

  void push_back(const T& element) {
    if (deq_.empty()) {
      T* arr;
      deq_.push_back(arr);
      deq_[0] = alloc_traits::allocate(alloc_, kColumn);
      alloc_traits::construct(alloc_, deq_[0], element);
      ++size_;
      return;
    }
    if (end_ind_2_ == kColumn - 1 && end_ind_1_ == deq_.size() - 1) {
      try {
        T* arr;
        deq_.push_back(arr);
        ++end_ind_1_;
        end_ind_2_ = 0;
        deq_[end_ind_1_] = alloc_traits::allocate(alloc_, kColumn);
        alloc_traits::construct(alloc_, deq_[end_ind_1_] + end_ind_2_, element);
        ++size_;
        return;
      } catch (...) {
        alloc_traits::destroy(alloc_, deq_[end_ind_1_] + end_ind_2_);
        alloc_traits::deallocate(alloc_, deq_[end_ind_1_], kColumn);
        --end_ind_1_;
        end_ind_2_ = kColumn - 1;
        throw;
      }
    }
    if (end_ind_2_ == kColumn - 1) {
      ++end_ind_1_;
      end_ind_2_ = 0;
    } else {
      ++end_ind_2_;
    }
    try {
      alloc_traits::construct(alloc_, deq_[end_ind_1_] + end_ind_2_, element);
    } catch (...) {
      alloc_traits::destroy(alloc_, deq_[end_ind_1_] + end_ind_2_);
      throw;
    }
    ++size_;
  }

  void push_back(T&& element) {
    if (deq_.empty()) {
      T* arr;
      deq_.push_back(arr);
      deq_[0] = alloc_traits::allocate(alloc_, kColumn);
      alloc_traits::construct(alloc_, deq_[0], std::move(element));
      ++size_;
      return;
    }
    if (end_ind_2_ == kColumn - 1 && end_ind_1_ == deq_.size() - 1) {
      try {
        T* arr;
        deq_.push_back(arr);
        ++end_ind_1_;
        end_ind_2_ = 0;
        deq_[end_ind_1_] = alloc_traits::allocate(alloc_, kColumn);
        alloc_traits::construct(alloc_, deq_[end_ind_1_] + end_ind_2_,
                                std::move(element));
        ++size_;
        return;
      } catch (...) {
        alloc_traits::destroy(alloc_, deq_[end_ind_1_] + end_ind_2_);
        alloc_traits::deallocate(alloc_, deq_[end_ind_1_], kColumn);
        --end_ind_1_;
        end_ind_2_ = kColumn - 1;
        throw;
      }
    }
    if (end_ind_2_ == kColumn - 1) {
      ++end_ind_1_;
      end_ind_2_ = 0;
    } else {
      ++end_ind_2_;
    }
    try {
      alloc_traits::construct(alloc_, deq_[end_ind_1_] + end_ind_2_,
                              std::move(element));
    } catch (...) {
      alloc_traits::destroy(alloc_, deq_[end_ind_1_] + end_ind_2_);
      throw;
    }
    ++size_;
  }

  void pop_back() {
    alloc_traits::destroy(alloc_, deq_[end_ind_1_] + end_ind_2_);
    --size_;
    if (end_ind_2_ == 0) {
      --end_ind_1_;
      end_ind_2_ = kColumn;
    }
    --end_ind_2_;
  }

  void pop_front() {
    alloc_traits::destroy(alloc_, deq_[bedin_ind_1_] + begin_ind_2_);
    --size_;
    if (begin_ind_2_ == kColumn - 1) {
      begin_ind_2_ = 0;
      ++bedin_ind_1_;
    } else {
      ++begin_ind_2_;
    }
  }

  void push_front(const T& element) {
    if (begin_ind_2_ == 0 && bedin_ind_1_ == 0) {
      try {
        T* arr;
        deq_.insert(deq_.begin(), arr);
        deq_[0] = alloc_traits::allocate(alloc_, kColumn);
        alloc_traits::construct(alloc_, deq_[0] + (kColumn - 1), element);
      } catch (...) {
        alloc_traits::destroy(alloc_, deq_[0] + (kColumn - 1));
        alloc_traits::deallocate(alloc_, deq_[0], kColumn);
        throw;
      }
      ++end_ind_1_;
      begin_ind_2_ = kColumn - 1;
      ++size_;
      return;
    }
    if (begin_ind_2_ == 0 && bedin_ind_1_ > 0) {
      --bedin_ind_1_;
      begin_ind_2_ = kColumn;
    }
    try {
      --begin_ind_2_;
      alloc_traits::construct(alloc_, deq_[bedin_ind_1_] + begin_ind_2_,
                              element);
    } catch (...) {
      alloc_traits::destroy(alloc_, deq_[bedin_ind_1_] + begin_ind_2_);
      ++begin_ind_2_;
      throw;
    }
    ++size_;
  }

  void push_front(T&& element) {
    if (begin_ind_2_ == 0 && bedin_ind_1_ == 0) {
      try {
        T* arr;
        deq_.insert(deq_.begin(), arr);
        deq_[0] = alloc_traits::allocate(alloc_, kColumn);
        alloc_traits::construct(alloc_, deq_[0] + (kColumn - 1),
                                std::move(element));
      } catch (...) {
        alloc_traits::destroy(alloc_, deq_[0] + (kColumn - 1));
        alloc_traits::deallocate(alloc_, deq_[0], kColumn);
        throw;
      }
      ++end_ind_1_;
      begin_ind_2_ = kColumn - 1;
      ++size_;
      return;
    }
    if (begin_ind_2_ == 0 && bedin_ind_1_ > 0) {
      --bedin_ind_1_;
      begin_ind_2_ = kColumn;
    }
    try {
      --begin_ind_2_;
      alloc_traits::construct(alloc_, deq_[bedin_ind_1_] + begin_ind_2_,
                              std::move(element));
    } catch (...) {
      alloc_traits::destroy(alloc_, deq_[bedin_ind_1_] + begin_ind_2_);
      ++begin_ind_2_;
      throw;
    }
    ++size_;
  }

  template <typename... Args>
  void emplace_back(Args&&... args) {
    if (deq_.empty()) {
      T* arr = alloc_traits::allocate(alloc_, kColumn);
      deq_.push_back(arr);
      alloc_traits::construct(alloc_, deq_[0], std::forward<Args>(args)...);
      ++size_;
      return;
    }
    if (end_ind_2_ == kColumn - 1 && end_ind_1_ == deq_.size() - 1) {
      try {
        T* arr = alloc_traits::allocate(alloc_, kColumn);
        deq_.push_back(arr);
        ++end_ind_1_;
        end_ind_2_ = 0;
        alloc_traits::construct(alloc_, deq_[end_ind_1_] + end_ind_2_,
                                std::forward<Args>(args)...);
        ++size_;
        return;
      } catch (...) {
        alloc_traits::destroy(alloc_, deq_[end_ind_1_] + end_ind_2_);
        alloc_traits::deallocate(alloc_, deq_[end_ind_1_], kColumn);
        --end_ind_1_;
        end_ind_2_ = kColumn - 1;
        throw;
      }
    }
    if (end_ind_2_ == kColumn - 1) {
      ++end_ind_1_;
      end_ind_2_ = 0;
    } else {
      ++end_ind_2_;
    }
    try {
      alloc_traits::construct(alloc_, deq_[end_ind_1_] + end_ind_2_,
                              std::forward<Args>(args)...);
    } catch (...) {
      alloc_traits::destroy(alloc_, deq_[end_ind_1_] + end_ind_2_);
      throw;
    }
    ++size_;
  }

  template <typename... Args>
  void emplace_front(Args&&... args) {
    if (begin_ind_2_ == 0 && bedin_ind_1_ == 0) {
      try {
        T* arr = alloc_traits::allocate(alloc_, kColumn);
        deq_.insert(deq_.begin(), arr);
        alloc_traits::construct(alloc_, deq_[0] + (kColumn - 1),
                                std::forward<Args>(args)...);
      } catch (...) {
        alloc_traits::destroy(alloc_, deq_[0] + (kColumn - 1));
        alloc_traits::deallocate(alloc_, deq_[0], kColumn);
        throw;
      }
      ++end_ind_1_;
      begin_ind_2_ = kColumn - 1;
      ++size_;
      return;
    }
    if (begin_ind_2_ == 0 && bedin_ind_1_ > 0) {
      --bedin_ind_1_;
      begin_ind_2_ = kColumn;
    }
    try {
      --begin_ind_2_;
      alloc_traits::construct(alloc_, deq_[bedin_ind_1_] + begin_ind_2_,
                              std::forward<Args>(args)...);
    } catch (...) {
      alloc_traits::destroy(alloc_, deq_[bedin_ind_1_] + begin_ind_2_);
      ++begin_ind_2_;
      throw;
    }
    ++size_;
  }

  template <bool IsConst>
  class CommonIterator {
   public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::conditional_t<IsConst, const T, T>;
    using value_type_vector =
        std::conditional_t<IsConst, const std::vector<T*>, std::vector<T*>>;
    using pointer = value_type*;
    using reference = value_type&;
    using difference_type = std::ptrdiff_t;

    CommonIterator(value_type_vector* vec, int nnn, int mmm) {
      vec_ = vec;
      idx_in_array_ = mmm;
      array_ = nnn;
    }

    reference operator*() const { return (*vec_)[array_][idx_in_array_]; }
    pointer operator->() const { return &((*vec_)[array_][idx_in_array_]); }

    CommonIterator& operator++() {
      if (idx_in_array_ == kColumnIter - 1) {
        idx_in_array_ = 0;
        ++array_;
        return *this;
      }
      ++idx_in_array_;
      return *this;
    }

    CommonIterator operator++(int) {
      CommonIterator<IsConst> copy = *this;
      ++(*this);
      return (copy);
    }

    CommonIterator& operator--() {
      if (idx_in_array_ == 0) {
        idx_in_array_ = kColumnIter - 1;
        --array_;
        return *this;
      }
      --idx_in_array_;
      return *this;
    }

    CommonIterator operator--(int) {
      CommonIterator<IsConst> copy = *this;
      --(*this);
      return (copy);
    }

    CommonIterator operator+(int n) {
      CommonIterator<IsConst> tmp = *this;
      if (n < 0) {
        return tmp - (-n);
      }
      tmp.idx_in_array_ += n;
      if (tmp.idx_in_array_ >= kColumnIter) {
        ++tmp.array_;
        tmp.idx_in_array_ %= kColumnIter;
      }
      return tmp;
    }

    CommonIterator& operator+=(int n) {
      *this = *this + n;
      return *this;
    }

    CommonIterator operator-(int n) {
      CommonIterator<IsConst> tmp = *this;
      if (n < 0) {
        return tmp + (-n);
      }
      if (static_cast<int>(tmp.idx_in_array_) < n) {
        n -= tmp.idx_in_array_;
        // tmp.mesto_v_massive_ = 0;
        tmp.idx_in_array_ = (kColumnIter * n - n) % kColumnIter;
        tmp.array_ -= (1 + (n - 1) / kColumnIter);
        return tmp;
      }
      tmp.idx_in_array_ -= n;
      return tmp;
    }

    CommonIterator& operator-=(int n) { return (*this = (*this - n)); }

    difference_type operator-(const CommonIterator& other) const {
      return (array_ - other.array_) * kColumnIter +
             (idx_in_array_ - other.idx_in_array_);
    }

    bool operator<(const CommonIterator& other) const {
      if (array_ == other.array_) {
        return idx_in_array_ < other.idx_in_array_;
      }
      return array_ < other.array_;
    }

    bool operator>(const CommonIterator& other) const { return other < *this; }

    bool operator!=(const CommonIterator& other) const {
      return !(*this == other);
    }

    bool operator==(const CommonIterator& other) const {
      return !(*this < other) && !(other < *this);
    }

    bool operator<=(const CommonIterator& other) const {
      return !(*this > other);
    }

    bool operator>=(const CommonIterator& other) const {
      return !(*this < other);
    }

   private:
    value_type_vector* vec_;
    static const int kColumnIter = 10000;
    int idx_in_array_;
    int array_;
  };

  using iterator = CommonIterator<false>;
  using const_iterator = CommonIterator<true>;

  template <typename Iterator>
  struct Reverseiterator {
   public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = typename Iterator::value_type;
    using pointer = typename Iterator::pointer;
    using reference = typename Iterator::reference;
    using difference_type = typename Iterator::difference_type;

    Reverseiterator(const Iterator& iter) : iter_(iter) {}
    Reverseiterator(const Reverseiterator& other) = default;
    Reverseiterator& operator=(const Reverseiterator& other) = default;
    reference operator*() const { return *iter_; }
    pointer operator->() const { return &(*iter_); }

    Reverseiterator<Iterator>& operator--() {
      ++iter_;
      return *this;
    }

    Reverseiterator operator--(int) {
      Reverseiterator copy = *this;
      --(*this);
      return copy;
    }

    Reverseiterator operator++(int) {
      Reverseiterator copy = *this;
      ++(*this);
      return copy;
    }

    Reverseiterator<Iterator>& operator++() {
      --iter_;
      return *this;
    }

    Reverseiterator operator+(int n) const {
      Iterator copy = iter_;
      copy -= n;
      return Reverseiterator(copy);
    }

    Reverseiterator& operator+=(int n) { return (*this = (*this + n)); }

    Reverseiterator operator-(int n) const {
      Iterator copy = iter_;
      copy += n;
      return Reverseiterator(copy);
    }

    Reverseiterator& operator-=(int n) { return (*this = (*this - n)); }

    difference_type operator-(const Reverseiterator& other) const {
      return (other.iter_ - iter_);
    }

    bool operator==(const Reverseiterator& other) const {
      return iter_ == other.iter_;
    }

    bool operator!=(const Reverseiterator& other) const {
      return iter_ != other.iter_;
    }

    bool operator<(const Reverseiterator& other) const {
      return iter_ > other.iter_;
    }

    bool operator>(const Reverseiterator& other) const {
      return iter_ < other.iter_;
    }

    bool operator<=(const Reverseiterator& other) const {
      return !(*this > other);
    }

    bool operator>=(const Reverseiterator& other) const {
      return !(*this < other);
    }

   private:
    Iterator iter_;
  };

  using reverse_iterator = Reverseiterator<iterator>;
  using const_reverse_iterator = Reverseiterator<const_iterator>;

  const_iterator begin() const {
    return const_iterator(&deq_, bedin_ind_1_, begin_ind_2_);
  }

  iterator begin() { return iterator(&deq_, bedin_ind_1_, begin_ind_2_); }

  const_iterator end() const {
    if (deq_.empty()) {
      return begin();
    }
    return const_iterator(
        &deq_, end_ind_1_ + static_cast<int>(end_ind_2_ == kColumn - 1),
        (end_ind_2_ + 1) % kColumn);
  }

  iterator end() {
    if (deq_.empty()) {
      return begin();
    }
    return iterator(&deq_,
                    end_ind_1_ + static_cast<int>(end_ind_2_ == kColumn - 1),
                    (end_ind_2_ + 1) % kColumn);
  }

  const_iterator cbegin() const {
    return const_iterator(&deq_, bedin_ind_1_, begin_ind_2_);
  }

  const_iterator cend() const {
    if (deq_.empty()) {
      return cbegin();
    }
    return const_iterator(
        &deq_, end_ind_1_ + static_cast<int>(end_ind_2_ == kColumn - 1),
        (end_ind_2_ + 1) % kColumn);
  }

  reverse_iterator rbegin() { return reverse_iterator(end() - 1); }
  reverse_iterator rbegin() const { return reverse_iterator(end() - 1); }
  reverse_iterator rend() { return reverse_iterator(begin() - 1); }
  reverse_iterator rend() const { return reverse_iterator(begin() - 1); }

  const_reverse_iterator crbegin() const {
    return const_reverse_iterator(cend() - 1);
  }

  const_reverse_iterator crend() const {
    return const_reverse_iterator(cbegin() - 1);
  }

  void insert(iterator itt, const T& element) {
    T temp = element;
    while (itt < end()) {
      T tmp = *itt;
      *itt = temp;
      temp = tmp;
      ++itt;
    }
    if (itt == end()) {
      push_back(temp);
    }
  }

  void emplace(iterator itt, T&& element) {
    T temp = std::move(element);
    while (itt < end()) {
      T tmp = *itt;
      *itt = temp;
      temp = tmp;
      ++itt;
    }
    if (itt == end()) {
      push_back(temp);
    }
  }

  void erase(iterator itt) {
    if (itt == begin()) {
      pop_front();
      return;
    }
    while (itt < end() - 1) {
      T tmp = *(itt + 1);
      *itt = tmp;
      ++itt;
    }
    if (itt == end() - 1) {
      pop_back();
    }
  }

  Allocator get_allocator() const { return alloc_; }

 private:
  std::vector<T*> deq_;
  const size_t kColumn = 10000;
  size_t size_ = 0;
  size_t end_ind_1_ = 0;
  size_t end_ind_2_ = 0;
  size_t bedin_ind_1_ = 0;
  size_t begin_ind_2_ = 0;
  Allocator alloc_ = Allocator();
};
