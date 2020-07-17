class RegionNotEvenException(Exception):
    def __init__(self, length):
        super().__init__()
        self.length = length

    def __repr__(self):
        return "RegionNotEvenException: The region length is not even, so the +1 nt position cannot be determined. Length: " + self.length

def verify_region_length_is_even(region_length):
    # If the region is an odd number, exit
    if region_length % 2 != 0:
        raise RegionNotEvenException(region_length)