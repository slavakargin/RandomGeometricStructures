from rgs.ribbons.Ribbon import (Ribbon, squares2ribbon)

def test_Ribbon():
    ribbon = Ribbon(0, 0, [0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1]) 
    assert ribbon.__str__() == "(0,0), [0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1]"


def test_squares2ribbon():
    squares = [(0, 2), (0, 3), (0, 4)]
    ribbon = squares2ribbon(squares)
    assert ribbon.__str__() == "(0,2), [1, 1]"
    assert ribbon.contains(0,4)
    
test_Ribbon()
test_squares2ribbon()