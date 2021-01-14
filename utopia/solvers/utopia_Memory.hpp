#ifndef UTOPIA_MEMORY_BASE_HPP
#define UTOPIA_MEMORY_BASE_HPP

#include "utopia_Layout.hpp"
#include "utopia_SolutionStatus.hpp"

namespace utopia {

template <class Vector>
class MemoryInterface {
 public:
  using SizeType = typename Traits<Vector>::SizeType;
  using Communicator = typename Traits<Vector>::Communicator;
  using Layout = typename Traits<Vector>::Layout;

  virtual ~MemoryInterface() = default;

  virtual void init_memory(const Layout &l) = 0;
};
}  // namespace utopia

#endif  // UTOPIA_MEMORY_BASE_HPP
