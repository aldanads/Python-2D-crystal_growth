# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 10:03:02 2023

@author: ALDANADS
"""

class Node:
    def __init__(self, data):
        self.data = data
        self.left = None
        self.right = None

def build_tree(arr):
    if not arr:
        return None
    if len(arr) == 1:
        return Node(arr[0])
    mid = len(arr)//2
    root = Node(None)
    #print(len(arr)/2)
    #print(arr[:mid], arr[mid:])
    root.left = build_tree(arr[:mid])
    root.right = build_tree(arr[mid:])
    return root

def update_data(root):
    if root is None:
        return
    
    # We start at the leafs
    update_data(root.left)
    update_data(root.right)
    
    # Superior nodes are the sum of their childrens - Leaf are tuples, but the
    # nodes are floats
    if root.left is not None and root.right is not None:
        if type(root.left.data) != tuple: 
            aux_l = root.left.data
        else:
            aux_l = root.left.data[0]
            
        if type(root.right.data) != tuple:
            aux_r = root.right.data
        else:
            aux_r = root.right.data[0]

        root.data = aux_l + aux_r

    elif root.left is not None:
        if type(root.left.data) != tuple: 
            root.data = root.left.data
        else: 
            root.data = root.left.data[0] 
            
    elif root.right is not None:
        if type(root.right.data) != tuple:
            root.data = root.right.data
        else: 
            root.data = root.right.data[0] 
    return root.data

def search_value(root,value,found = False):
    if root is None or found is not False:
        return found

    if (root.left is None) and (root.right is None) and (root.data == value):
        found = root.data
    
    found = search_value(root.left, value,found)
    found = search_value(root.right, value,found)
    
    return found

# example usage:
# arr = [1, 2, 3, 4, 5, 6, 7]
# arr.sort(reverse = True)
# root = build_tree(arr)
# total = update_data(root)
# print(total)