function print_title(title::String)
    println(title)
    for _ in 1:length(title) print("-") end
    print("\n\n")
    return
end

function print_segment(symbol::String)
    print("#")
    for _ in 1:80
        print("$symbol")
    end
    print("\n")
    return
end

function print_text(text::String)
    len = length(text)
    remaining_space = 80 - len
    print("#")
    for _ in 1:fld(remaining_space, 2)
        print(" ")
    end
    print(text)
    for _ in 1:(fld(remaining_space, 2) + 1)
        print(" ")
    end
    print("\n")
    return
end

function print_header()
    print_segment("-")
    print_segment(" ")
    print_text("Bilevel DNEP")
    print_segment(" ")
    print_text("Graduation work")
    print_segment(" ")
    print_segment("-")
    println("# @ Manon Cornet\n")
    return
end