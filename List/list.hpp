#pragma once
#include <iostream>
#include <memory>

template <typename T, typename Allocator = std::allocator<T>>
class List {
 private:
  struct BaseNode {
    BaseNode* prev;
    BaseNode* next;
    BaseNode() = default;
    BaseNode(BaseNode* prev, BaseNode* next) : prev(prev), next(next) {}
  };

  struct Node : BaseNode {
    T value;
    Node() = default;
    Node(const T& value) : value(value) {}
    Node(BaseNode* prev, BaseNode* next, T& value)
        : BaseNode(prev, next), value(value) {}
    Node(T&& value) : value(std::move(value)) {}
  };

 public:
  using value_type = std::remove_cv_t<T>;
  using allocator_type = Allocator;
  using alloc_traits = std::allocator_traits<Allocator>;
  using node_alloc = typename alloc_traits::template rebind_alloc<Node>;
  using node_alloc_traits = typename alloc_traits::template rebind_traits<Node>;

 private:
  BaseNode fakenode_;
  size_t size_ = 0;
  node_alloc alloc_;

 public:
  List() {}
  List(size_t count, const T& value, const Allocator& alloc = Allocator())
      : alloc_(alloc) {
    Node* tmp;
    size_t iii = 0;
    try {
      for (; iii < count; ++iii) {
        tmp = node_alloc_traits::allocate(alloc_, 1);
        node_alloc_traits::construct(alloc_, tmp, value);
        if (iii == 0) {
          tmp->prev = &fakenode_;
          tmp->next = &fakenode_;
          fakenode_.next = static_cast<BaseNode*>(tmp);
          fakenode_.prev = static_cast<BaseNode*>(tmp);
        } else {
          fakenode_.prev->next = static_cast<BaseNode*>(tmp);
          tmp->prev = fakenode_.prev;
          fakenode_.prev = static_cast<BaseNode*>(tmp);
          tmp->next = &fakenode_;
        }
      }
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      for (size_t j = 0; j < iii; ++j) {
        pop_back();
      }
      throw;
    }
    size_ = count;
  }

  explicit List(size_t count, const Allocator& alloc = Allocator())
      : alloc_(alloc) {
    Node* tmp;
    size_t iii = 0;
    try {
      for (; iii < count; ++iii) {
        tmp = node_alloc_traits::allocate(alloc_, 1);
        node_alloc_traits::construct(alloc_, tmp);
        if (iii == 0) {
          tmp->prev = &fakenode_;
          tmp->next = &fakenode_;
          fakenode_.next = static_cast<BaseNode*>(tmp);
          fakenode_.prev = static_cast<BaseNode*>(tmp);
        } else {
          fakenode_.prev->next = static_cast<BaseNode*>(tmp);
          tmp->prev = fakenode_.prev;
          fakenode_.prev = static_cast<BaseNode*>(tmp);
          tmp->next = &fakenode_;
        }
      }
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      for (size_t j = 0; j < iii; ++j) {
        pop_back();
      }
      throw;
    }
    size_ = count;
  }

  List(const List& other) {
    alloc_ =
        node_alloc_traits::select_on_container_copy_construction(other.alloc_);
    auto itt = other.begin();
    Node* tmp;
    try {
      for (; itt != other.end(); ++itt) {
        tmp = node_alloc_traits::allocate(alloc_, 1);
        node_alloc_traits::construct(alloc_, tmp, *itt);
        if (itt == other.begin()) {
          tmp->prev = &fakenode_;
          tmp->next = &fakenode_;
          fakenode_.next = tmp;
          fakenode_.prev = tmp;
        } else {
          fakenode_.prev->next = tmp;
          tmp->prev = fakenode_.prev;
          fakenode_.prev = tmp;
          tmp->next = &fakenode_;
        }
      }
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      for (auto itt2 = other.begin(); itt2 != itt; ++itt2) {
        pop_back();
      }
      throw;
    }
    size_ = other.size_;
  }

  List(const List& other, const Allocator& alloc) {
    auto itt = other.begin();
    Node* tmp;
    try {
      for (; itt != other.end(); ++itt) {
        tmp = node_alloc_traits::allocate(alloc, 1);
        node_alloc_traits::construct(alloc, tmp, *itt);
        if (other.size_ == 1) {
          tmp->prev = &fakenode_;
          tmp->next = &fakenode_;
          fakenode_.next = tmp;
          fakenode_.prev = tmp;
        } else {
          fakenode_.prev->next = tmp;
          tmp->prev = fakenode_.prev;
          fakenode_.prev = tmp;
          tmp->next = &fakenode_;
        }
      }
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      for (auto itt2 = other.begin(); itt2 != itt; ++itt2) {
        pop_back();
      }
      throw;
    }
    alloc_ = alloc;
    size_ = other.size_;
  }

  List(std::initializer_list<T> init, const Allocator& alloc = Allocator()) {
    alloc_ = alloc;
    auto itt = init.begin();
    Node* tmp;
    try {
      for (; itt != init.end(); ++itt) {
        tmp = node_alloc_traits::allocate(alloc_, 1);
        node_alloc_traits::construct(alloc_, tmp, *itt);
        if (itt == init.begin()) {
          tmp->prev = &fakenode_;
          tmp->next = &fakenode_;
          fakenode_.next = tmp;
          fakenode_.prev = tmp;
        } else {
          fakenode_.prev->next = tmp;
          tmp->prev = fakenode_.prev;
          fakenode_.prev = tmp;
          tmp->next = &fakenode_;
        }
      }
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      for (auto itt2 = init.begin(); itt2 != itt; ++itt2) {
        pop_back();
      }
      throw;
    }
    size_ = init.size();
  }

  ~List() {
    for (size_t i = 0; i < size_; ++i) {
      Node* tmp = static_cast<Node*>(fakenode_.prev->prev);
      node_alloc_traits::destroy(alloc_, static_cast<Node*>(fakenode_.prev));
      node_alloc_traits::deallocate(alloc_, static_cast<Node*>(fakenode_.prev),
                                    1);
      tmp->next = &fakenode_;
      fakenode_.prev = tmp;
    }
  }
  List& operator=(const List& other) {
    List tmp(other);
    std::swap(fakenode_, tmp.fakenode_);
    std::swap(size_, tmp.size_);
    if (std::allocator_traits<decltype(
            other.alloc_)>::propagate_on_container_copy_assignment::value) {
      alloc_ = other.alloc_;
      return *this;
    }
    return *this;
  }

  T& front() { return static_cast<Node*>(fakenode_.next)->value; }
  const T& front() const { return static_cast<Node*>(fakenode_.next)->value; }
  T& back() { return static_cast<Node*>(fakenode_.prev)->value; }
  const T& back() const { return static_cast<Node*>(fakenode_.prev)->value; }
  bool empty() const { return size_ == 0; }
  size_t size() const { return size_; }

  void push_back(const T& elem) {
    Node* tmp;
    try {
      tmp = node_alloc_traits::allocate(alloc_, 1);
      node_alloc_traits::construct(alloc_, tmp, elem);
      if (empty()) {
        fakenode_.prev = static_cast<BaseNode*>(tmp);
        fakenode_.next = static_cast<BaseNode*>(tmp);
        tmp->prev = &fakenode_;
        tmp->next = &fakenode_;
      } else {
        tmp->prev = fakenode_.prev;
        tmp->next = &fakenode_;
        fakenode_.prev->next = static_cast<BaseNode*>(tmp);
        fakenode_.prev = static_cast<BaseNode*>(tmp);
      }
      ++size_;
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      throw;
    }
  }

  void push_back(T&& elem) {
    Node* tmp;
    try {
      tmp = node_alloc_traits::allocate(alloc_, 1);
      node_alloc_traits::construct(alloc_, tmp, std::move(elem));
      if (empty()) {
        fakenode_.prev = static_cast<BaseNode*>(tmp);
        fakenode_.next = static_cast<BaseNode*>(tmp);
        tmp->prev = &fakenode_;
        tmp->next = &fakenode_;
      } else {
        tmp->prev = fakenode_.prev;
        tmp->next = &fakenode_;
        fakenode_.prev->next = static_cast<BaseNode*>(tmp);
        fakenode_.prev = static_cast<BaseNode*>(tmp);
      }
      ++size_;
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      throw;
    }
  }

  void pop_back() {
    Node* tmp = static_cast<Node*>(fakenode_.prev->prev);
    node_alloc_traits::destroy(alloc_, static_cast<Node*>(fakenode_.prev));
    node_alloc_traits::deallocate(alloc_, static_cast<Node*>(fakenode_.prev),
                                  1);
    tmp->next = &fakenode_;
    fakenode_.prev = static_cast<BaseNode*>(tmp);
    --size_;
  }

  void push_front(const T& elem) {
    Node* tmp;
    try {
      tmp = node_alloc_traits::allocate(alloc_, 1);
      node_alloc_traits::construct(alloc_, tmp, elem);
      if (empty()) {
        fakenode_.prev = static_cast<BaseNode*>(tmp);
        fakenode_.next = static_cast<BaseNode*>(tmp);
        tmp->prev = &fakenode_;
        tmp->next = &fakenode_;
      } else {
        tmp->prev = &fakenode_;
        tmp->next = fakenode_.next;
        fakenode_.next->prev = static_cast<BaseNode*>(tmp);
        fakenode_.next = static_cast<BaseNode*>(tmp);
      }
      ++size_;
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      throw;
    }
  }

  void push_front(T&& elem) {
    Node* tmp;
    try {
      tmp = node_alloc_traits::allocate(alloc_, 1);
      node_alloc_traits::construct(alloc_, tmp, std::move(elem));
      if (empty()) {
        fakenode_.prev = static_cast<BaseNode*>(tmp);
        fakenode_.next = static_cast<BaseNode*>(tmp);
        tmp->prev = &fakenode_;
        tmp->next = &fakenode_;
      } else {
        tmp->prev = &fakenode_;
        tmp->next = fakenode_.next;
        fakenode_.next->prev = static_cast<BaseNode*>(tmp);
        fakenode_.next = static_cast<BaseNode*>(tmp);
      }
      ++size_;
    } catch (...) {
      node_alloc_traits::deallocate(alloc_, tmp, 1);
      throw;
    }
  }

  void pop_front() {
    Node* tmp = static_cast<Node*>(fakenode_.next->next);
    node_alloc_traits::destroy(alloc_, static_cast<Node*>(fakenode_.next));
    node_alloc_traits::deallocate(alloc_, static_cast<Node*>(fakenode_.next),
                                  1);
    tmp->prev = &fakenode_;
    fakenode_.next = static_cast<BaseNode*>(tmp);
    --size_;
  }

  template <bool IsConst>
  class CommonIterator {
   private:
    BaseNode* node_;

   public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = std::conditional_t<IsConst, const T, T>;
    using node_type = std::conditional_t<IsConst, const BaseNode, BaseNode>;
    using pointer = value_type*;
    using reference = value_type&;
    using difference_type = std::ptrdiff_t;

    CommonIterator(node_type* node) { node_ = node; }
    reference operator*() const { return static_cast<Node*>(node_)->value; }
    pointer operator->() const { return &(static_cast<Node*>(node_)->value); }

    CommonIterator& operator++() {
      node_ = node_->next;
      return *this;
    }

    CommonIterator operator++(int) {
      CommonIterator copy = *this;
      ++(*this);
      return (copy);
    }

    CommonIterator& operator--() {
      node_ = node_->prev;
      return *this;
    }

    CommonIterator operator--(int) {
      CommonIterator copy = *this;
      --(*this);
      return (copy);
    }

    CommonIterator operator+(int n) {
      CommonIterator tmp = *this;
      if (n < 0) {
        return tmp - (-n);
      }
      for (int i = 0; i < n; ++i) {
        ++tmp.node_;
      }
      return tmp;
    }

    CommonIterator& operator+=(int n) {
      *this = *this + n;
      return *this;
    }

    CommonIterator operator-(int n) {
      CommonIterator tmp = *this;
      if (n < 0) {
        return tmp + (-n);
      }
      for (int i = 0; i < n; ++i) {
        --tmp;
      }
      return tmp;
    }

    CommonIterator& operator-=(int n) { return (*this = (*this - n)); }

    bool operator!=(const CommonIterator& other) const {
      return node_ != other.node_;
    }

    bool operator==(const CommonIterator& other) const {
      return node_ == other.node_;
    }
  };

  using iterator = CommonIterator<false>;
  using const_iterator = CommonIterator<true>;

  template <typename Iterator>
  struct Reverseiterator {
   public:
    using iterator_category = std::bidirectional_iterator_tag;
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

    bool operator==(const Reverseiterator& other) const {
      return iter_ == other.iter_;
    }

    bool operator!=(const Reverseiterator& other) const {
      return iter_ != other.iter_;
    }

   private:
    Iterator iter_;
  };

  using reverse_iterator = Reverseiterator<iterator>;
  using const_reverse_iterator = Reverseiterator<const_iterator>;

  iterator begin() const { return iterator(fakenode_.next); }
  iterator begin() { return iterator(fakenode_.next); }
  iterator end() const { return iterator(const_cast<BaseNode*>(&fakenode_)); }
  iterator end() { return iterator(&fakenode_); }
  const_iterator cbegin() const { return const_iterator(fakenode_.next); }
  const_iterator cend() const { return const_iterator(&fakenode_); }
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

  Allocator get_allocator() const { return alloc_; }
};
