// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_bvh.h"
#include <numeric> // iota

// Adapted from:
// https://github.com/mattoverby/mclscene/blob/master/include/MCL/BVH.hpp

namespace admmpd {

template <typename T, int DIM>
void AABBTree<T,DIM>::clear()
{
    root = std::make_shared<Node>();
}

template <typename T, int DIM>
void AABBTree<T,DIM>::init(const std::vector<AABB> &leaves)
{
    root = std::make_shared<Node>();
    int np = leaves.size();
    if (np==0)
        return;
    std::vector<int> queue(np);
    std::iota(queue.begin(), queue.end(), 0);
    create_children(root.get(), queue, leaves);
}

template <typename T, int DIM>
void AABBTree<T,DIM>::update(const std::vector<AABB> &leaves)
{
    if (!root || (int)leaves.size()==0)
        return;
    update_children(root.get(), leaves);
}

template <typename T, int DIM>
bool AABBTree<T,DIM>::traverse(Traverser<T,DIM> &traverser) const
{
    if (!root)
        return false;
    return traverse_children(root.get(), traverser);
}

// If we are traversing with function pointers, we'll just
// wrap our own traverser that calls these functions to
// avoid duplicate traverse_children code.
template <typename T, int DIM>
class TraverserFromFunctionPtrs : public Traverser<T,DIM>
{
using typename Traverser<T,DIM>::AABB;
public:
    std::function<void(const AABB&, bool&, const AABB&, bool&, bool&)> t;
    std::function<bool(const AABB&, int)> s;
    void traverse(
		const AABB &left_aabb, bool &go_left,
		const AABB &right_aabb, bool &go_right,
		bool &go_left_first)
    {
        t(left_aabb, go_left, right_aabb, go_right, go_left_first);
    }
	bool stop_traversing(const AABB &aabb, int prim)
    {
        return s(aabb,prim);
    }
};

template <typename T, int DIM>
bool AABBTree<T,DIM>::traverse(
    std::function<void(const AABB&, bool&, const AABB&, bool&, bool&)> t,
    std::function<bool(const AABB&, int)> s) const
{
    if (!root)
        return false;
    TraverserFromFunctionPtrs<T,DIM> traverser;
    traverser.t = t;
    traverser.s = s;
    return traverse_children(root.get(), traverser);
}

template <typename T, int DIM>
void AABBTree<T,DIM>::create_children(
    Node *node,
    std::vector<int> &queue,
    const std::vector<AABB> &leaves)
{
	node->aabb.setEmpty();
    int n_queue = queue.size();
	if (n_queue == 1)
    {
	    node->prims.emplace_back(queue[0]);
		node->aabb = leaves[queue[0]];
		return;
	}

    for (int i=0; i<n_queue; ++i)
    {
        int q = queue[i];
        node->aabb.extend(leaves[q]);
    }

    struct SortByAxis
    {
        int axis;
        const std::vector<AABB> &aabbs;
        SortByAxis(int axis_, const std::vector<AABB> &aabbs_) :
            axis(axis_), aabbs(aabbs_) {}
        bool operator()(size_t l, size_t r) const
        {
            return aabbs[l].center()[axis] < aabbs[r].center()[axis];
        }
    };

    // Sort tree and split queue
    int sort_axis = 0;
    VecType sizes = node->aabb.sizes();
    sizes.maxCoeff(&sort_axis);
    std::sort(queue.begin(), queue.end(), SortByAxis(sort_axis,leaves));
	std::vector<int> left_queue(queue.begin(), queue.begin()+(n_queue/2));
	std::vector<int> right_queue(queue.begin()+(n_queue/2), queue.end());

    // Recursive top-down constructrion
    node->left = new Node();
	create_children(node->left, left_queue, leaves);
    node->right = new Node();
	create_children(node->right, right_queue, leaves);

} // end create children

template <typename T, int DIM>
void AABBTree<T,DIM>::update_children(
    Node *node,
    const std::vector<AABB> &leaves)
{
	node->aabb.setEmpty();

    if (node->is_leaf())
    {
        int np = node->prims.size();
        for (int i=0; i<np; ++i)
        {
            node->aabb.extend(leaves[node->prims[i]]);
        }
        return;
    }

	if (node->left != nullptr)
    {
		update_children(node->left, leaves);
		node->aabb.extend(node->left->aabb);
	}
	if (node->right != nullptr)
    {
		update_children(node->right, leaves);
		node->aabb.extend(node->right->aabb);
	}

} // end update children

template <typename T, int DIM>
bool AABBTree<T,DIM>::traverse_children(
    const Node *node,
    Traverser<T,DIM> &traverser ) const
{
	if( node->is_leaf() ){
		int np = node->prims.size();
		for(int i=0; i<np; ++i)
        {
			if(traverser.stop_traversing(node->aabb, node->prims[i]))
                return true;
		}
		return false;
	}

	bool go_left = true;
	bool go_right = true;
	bool go_left_first = true;
    const AABB &left_aabb = (node->left == nullptr ? AABB() : node->left->aabb);
    const AABB &right_aabb = (node->right == nullptr ? AABB() : node->right->aabb);
	traverser.traverse(
		left_aabb, go_left,
		right_aabb, go_right,
		go_left_first );

	if (go_left && go_right)
    {
        if (go_left_first)
        {
            if (traverse_children(node->left, traverser)) { return true; }
		    else { return traverse_children(node->right, traverser); }
        }
        else
        {
            if (traverse_children(node->right, traverser)) { return true; }
		    else { return traverse_children(node->left, traverser); }
        }
	}
	if (go_left && !go_right)
    {
		return traverse_children(node->left, traverser);
	}
	if (!go_left && go_right)
    {
		return traverse_children(node->right, traverser);
	}

	return false;

} // end traverse children

// Compile types
template class admmpd::AABBTree<double,2>;
template class admmpd::AABBTree<double,3>;
template class admmpd::AABBTree<float,2>;
template class admmpd::AABBTree<float,3>;
template class admmpd::TraverserFromFunctionPtrs<double,2>;
template class admmpd::TraverserFromFunctionPtrs<double,3>;
template class admmpd::TraverserFromFunctionPtrs<float,2>;
template class admmpd::TraverserFromFunctionPtrs<float,3>;

} // namespace admmpd