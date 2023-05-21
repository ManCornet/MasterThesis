using PackageCompiler
import Pkg

all_packages = String[]
for (uuid, dep) in Pkg.dependencies()
    dep.is_direct_dep || continue
    dep.version === nothing && continue
    push!(all_packages,dep.name)
end

# Remove unneeded packages
do_not_include = ["PackageCompiler"]
package_list = filter(x -> x ∉ do_not_include, all_packages)

create_sysimage(package_list, sysimage_path="JuliaSysimage.so")