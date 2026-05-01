import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="TENGU: Transcript-signal ENrichment and Grouping Unit")
    subparsers = parser.add_subparsers(dest="command", help="TENGU subcommands")

    # (1) tengu segmentation
    seg_parser = subparsers.add_parser("segmentation", help="Aggregate 2µm bins into segmented units")
    seg_parser.add_argument("--counts", required=True, help="Path to counts folder")
    seg_parser.add_argument("--image", required=True, help="Path to tissue image")
    seg_parser.add_argument("--out_h5ad", required=True, help="Output h5ad file")
    seg_parser.add_argument("--out_cellmark", required=True, help="Output cellmark CSV")

    # (3) tengu annotation
    ann_parser = subparsers.add_parser("annotation", help="Annotate cell types")
    ann_parser.add_argument("--counts", required=True, help="Path to TENGU_counts")
    ann_parser.add_argument("--model", required=True, help="Path to CellTypist model")
    ann_parser.add_argument("--out", required=True, help="Output predicted cell type CSV")

    # (4) tengu communication
    com_parser = subparsers.add_parser("communication", help="Derive spatially aware cell-cell communication")
    com_parser.add_argument("--counts", required=True, help="Path to TENGU_counts")
    com_parser.add_argument("--database", required=True, help="Path to CellNEST database")
    com_parser.add_argument("--out", required=True, help="Output CCC results CSV")

    # (5) tengu simulation
    sim_parser = subparsers.add_parser("simulation", help="Probabilistic cell simulation framework")
    sim_parser.add_argument("--counts", required=True, help="Path to TENGU_counts")
    sim_parser.add_argument("--out_zarr", required=True, help="Output zarr object")
    sim_parser.add_argument("--out_cellmark", required=True, help="Output simulation cellmark CSV")

    args = parser.parse_args()

    if args.command == "segmentation":
        print(f"Running segmentation with: counts={args.counts}, image={args.image}")
    elif args.command == "annotation":
        print(f"Running annotation with: counts={args.counts}, model={args.model}")
    elif args.command == "communication":
        print(f"Running communication with: counts={args.counts}, database={args.database}")
    elif args.command == "simulation":
        print(f"Running simulation with: counts={args.counts}, out_zarr={args.out_zarr}")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
