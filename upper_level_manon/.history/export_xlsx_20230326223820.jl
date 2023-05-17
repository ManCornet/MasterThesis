
using OrderedCollections

function xlsx_output(model)
    
    dict_output = Dict(k => typeof(v) for (k, v) in object_dictionary(model) if v isa AbstractArray{VariableRef})

    # conditions for this preprocessing loop to work
    # - [K,L,...]  dimension K must be the 1st and dimension L the 2nd
    # - [K,...,Y,...] dimension K must be before dimension Y, in case both sets have the same length
    # - [...,Y,T,...] dimension T must follow Y
    alpha_name="alpha"
    for (key,value) in dict_output
        if ndims(value)>2
            value=value.data
            size_value=size(value)
            for (dim_idx,dim_length) in enumerate(size_value)
                if dim_length==length(K) && haskey(dict_output, alpha_name) && size_value[1:2]==size(dict_output[alpha_name])
                    dict_output[key]=dropdims(sum(value.*dict_output[alpha_name],dims=dim_idx),dims=dim_idx) # eliminates dimension K
                elseif dim_length==length(Y)
                    dict_output[key]=reshape(value, [size_value[1:dim_idx-1];size_value[dim_idx]*size_value[dim_idx+1]]) # eliminates dimension Y 
                end
            end
        end
    end


    XLSX.openxlsx("output "*(Dates.format(now(), "yyyy-mm-dd HH-MM-SS"))*".xlsx", mode="w") do xf
        sheet_n=1
        for (key,value) in dict_output
            sheet_n==1 && XLSX.rename!(xf[1], key)
            sheet_n>=2 && XLSX.addsheet!(xf, key)
            sheet = xf[sheet_n]
            # println(key, " ", value)
            if ndims(value)>=1
                items=axes(value)[1] # gets back the set L, Ns or Nu
                item="item"
                if issetequal(items,L)
                    item="line"
                elseif issetequal(items,Ns)
                    item="sub_node"
                elseif issetequal(items,Nu)
                    item="user_node"
                end
            end
            if ndims(value)==0
                df=DataFrame(Dict(key=>value))
                XLSX.writetable!(sheet, df)
            elseif ndims(value)==1
                df=DataFrame(Dict(key=>value.data))
                insertcols!(df,1,item=>1:size(value)[1]) # add 1st column: n° of node/line
                XLSX.writetable!(sheet, df)
            elseif ndims(value)==2
                items=axes(value)[1]
                steps=1:size(value)[2]            
                df=DataFrame(item=>items) # 1st column: n° of node/line
                [df[!,"step $i"]=value[:,i].data for i in steps] # add columns for each time step
                XLSX.writetable!(sheet, df)
            end
            sheet_n=sheet_n+1
        end
    end
end

# K L Y T