#ifndef _FUNCTIONAL_ARITHMETIC_HPP_
#define _FUNCTIONAL_ARITHMETIC_HPP_

#include <functional>

// Arithmetic function wrapper
template <class... Args>
class Function {
  std::function<double(Args&&...)> f;

public:
  explicit Function(const std::function<double(Args&&...)>& f) : f(f) {};

  auto operator()(Args&&... args) const {
    return f(std::forward<Args>(args)...);
  };
};


// Function addition
template <class... Args>
auto operator+(
  const Function<Args...>& expression1,
  const Function<Args...>& expression2) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression1(std::forward<Args>(args)...) + 
             expression2(std::forward<Args>(args)...);
    });
};

template <class... Args>
auto operator+(
  const Function<Args...>& expression,
  double scalar) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression(std::forward<Args>(args)...) + scalar;
    });
};

template <class... Args>
auto operator+(
  double scalar,
  const Function<Args&&...>& expression) {
  return expression + scalar;
};


// Function subtraction
template <class... Args>
auto operator-(
  const Function<Args...>& expression1,
  const Function<Args...>& expression2) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression1(std::forward<Args>(args)...) - 
             expression2(std::forward<Args>(args)...);
    });
};

template <class... Args>
auto operator-(
  const Function<Args...>& expression,
  double scalar) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression(std::forward<Args>(args)...) - scalar;
    });
};

template <class... Args>
auto operator-(
  double scalar,
  const Function<Args...>& expression) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return scalar - expression(std::forward<Args>(args)...);
    });
};


// Function multiplication
template <class... Args>
auto operator*(
  const Function<Args...>& expression1,
  const Function<Args...>& expression2) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression1(std::forward<Args>(args)...) * 
             expression2(std::forward<Args>(args)...);
    });
};

template <class... Args>
auto operator*(
  const Function<Args...>& expression,
  double scalar) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression(std::forward<Args>(args)...) * scalar;
    });
};

template <class... Args>
auto operator*(
  double scalar,
  const Function<Args...>& expression) {
  return expression * scalar;
};


// Function division
template <class... Args>
auto operator/(
  const Function<Args...>& expression1,
  const Function<Args...>& expression2) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression1(std::forward<Args>(args)...) / 
             expression2(std::forward<Args>(args)...);
    });
};

template <class... Args>
auto operator/(
  const Function<Args...>& expression,
  double scalar) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return expression(std::forward<Args>(args)...) / scalar;
    });
};

template <class... Args>
auto operator/(
  double scalar,
  const Function<Args...>& expression) {
  return Function<Args...>(
    [=] (Args&&... args) {
      return scalar / expression(std::forward<Args>(args)...);
    });
};

#endif // !_FUNCTIONAL_ARITHMETIC_HPP_