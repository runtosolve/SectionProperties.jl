struct MyCustomStruct
    a::Int
    b::Union{String, Nothing}

    # Inner constructor for when 'b' is provided
    MyCustomStruct(a::Int, b::String) = new(a, b)
    # Inner constructor for when 'b' is omitted (defaults to nothing)
    MyCustomStruct(a::Int) = new(a, nothing)
end

# Creating instances
instance6 = MyCustomStruct(10, "value")
instance7 = MyCustomStruct(20) # b defaults to nothing