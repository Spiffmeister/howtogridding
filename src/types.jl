
"""
    NodeType

- Used by `FirstDerivative()` and `SecondDerivative()` to indicate which node type should be used.
- `Left` and `Right` are used by the `SAT()` function to determine which side of the SAT to call.

Node types are:
- Left
- Internal
- Right
- Up
- Down
"""
struct NodeType{T,D} end
const Left = NodeType{:Left,1}()
const Internal = NodeType{:Internal,0}()
const Right = NodeType{:Right,1}()
const Up = NodeType{:Left,2}()
const Down = NodeType{:Right,2}()


struct DerivativeOrder{O} end
