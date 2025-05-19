# Load required libraries
library(geosphere)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)  # For statistical annotations

# 1. Load FST matrix (full dataset)
fst_data <- read.table(text = "Region NOR-KOF NOR-TRO Iceland NOR-CW NOR-CE NOR-SW NOR-SE Finland Lithuania Britain Germany Croatia Italy IBE-NE
NOR-ALT 0.77 0.69 0.76 0.58 0.53 0.62 0.52 0.40 0.45 0.62 0.52 0.54 0.61 0.74
NOR-KOF NA 0.60 0.75 0.52 0.45 0.58 0.43 0.41 0.41 0.60 0.45 0.46 0.57 0.73
NOR-TRO NA NA 0.28 0.17 0.25 0.20 0.30 0.39 0.40 0.26 0.28 0.35 0.28 0.32
Iceland NA NA NA 0.20 0.28 0.20 0.34 0.44 0.45 0.26 0.31 0.39 0.31 0.33
NOR-CW NA NA NA NA 0.04 0.03 0.09 0.22 0.24 0.08 0.09 0.18 0.12 0.21
NOR-CE NA NA NA NA NA 0.10 0.06 0.14 0.16 0.14 0.09 0.15 0.18 0.31
NOR-SW NA NA NA NA NA NA 0.14 0.28 0.29 0.11 0.14 0.23 0.16 0.26
NOR-SE NA NA NA NA NA NA NA 0.11 0.13 0.18 0.07 0.13 0.19 0.36
Finland NA NA NA NA NA NA NA NA 0.05 0.29 0.15 0.17 0.29 0.45
Lithuania NA NA NA NA NA NA NA NA NA 0.30 0.15 0.18 0.30 0.47
Britain NA NA NA NA NA NA NA NA NA NA 0.13 0.23 0.16 0.21
Germany NA NA NA NA NA NA NA NA NA NA NA 0.10 0.11 0.26
Croatia NA NA NA NA NA NA NA NA NA NA NA NA 0.18 0.36
Italy NA NA NA NA NA NA NA NA NA NA NA NA NA 0.19", 
                       header = TRUE, row.names = 1, check.names = FALSE)

# Define all 15 FST regions
fst_regions <- c("NOR-ALT", "NOR-KOF", "NOR-TRO", "Iceland", "NOR-CW", "NOR-CE", 
                 "NOR-SW", "NOR-SE", "Finland", "Lithuania", "Britain", "Germany", 
                 "Croatia", "Italy", "IBE-NE")
n <- length(fst_regions)  # 15

# 2. Load GPS coordinates (same as before, no changes needed)
gps_data <- read.table(text = "Latitude Longitude Region
47.0539432	13.0999286	Austria
44.770848	15.652355	Croatia
45.890571	15.963867	Croatia
45.888792	15.930374	Croatia
45.838877	15.598735	Croatia
45.687132	15.381353	Croatia
50.8685896	14.2971467	Czech_with_germany
55.5703	9.7466	Denmark
55.214605	11.4641154	Denmark
43.3099	−4.4697	IBE-NW
41.6176	0.62	IBE-NE
42.6712	0.9798	IBE-NE
43.5339	−6.5271	IBE-NW
43.1344	−4.888	IBE-NW
43.1425	−4.938	IBE-NW
43.0697	−5.0328	IBE-NW
41.5212	0.3531	IBE-NE
43.2504	−5.9833	IBE-NW
37.7796	−3.7849	IBE-S
40.2938	−5.0091	IBE-S
37.2614	−6.9447	IBE-S
37.891	−6.5609	IBE-S
42.5076	−0.6732	IBE-NE
42.8117	−1.6483	IBE-NE
42.6469	0.0349	IBE-NE
42.7383	−0.75	IBE-NE
42.6496	1.0347	IBE-NE
42.7557	−0.829	IBE-NE
60.4867	22.1534	Finland
60.7196	24.6824	Finland
60.3431	25.6795	Finland
60.4775	24.8444	Finland
60.1949	21.3349	Finland
60.5367	21.8766	Finland
60.2	21.4	Finland
60.15	21.59	Finland
60.2292	25.0233	Finland
60.7262	21.9093	Finland
60.4051	25.1562	Finland
60.4039	23.1252	Finland
60.627	25.8153	Finland
60.356	22.9651	Finland
60.1061	23.6782	Finland
60.0673	23.2879	Finland
60.2076	23.8066	Finland
60.4651	25.5254	Finland
61.8	29.3	Finland
59.8428	23.2446	Finland
60.3706	22.9796	Finland
60.0824	24.1592	Finland
42.0396042	9.0128926	France
48.801407	2.130122	France
50.015654	2.6973567	France
48.42353	7.66326	France
51.8236	11.2866	Germany
47.5557	10.0225	Germany
47.96666667	7.83333333	Germany
50.5833	8.65	Germany
50.6520515	9.1624376	Germany
50.5287	8.6842	Germany
51.1583	13.6814	Germany
51.4192	11.128	Germany
50.9847	11.3225	Germany
49.7853	11.3606	Germany
47.5857	10.5587	Germany
50.58	13	Germany
47.80331	8.0369	Germany
64.010806	−20.98075	Iceland
64.759279	−21.593767	Iceland
65.950069	−19.487877	Iceland
63.96063	−22.398191	Iceland
63.84166667	−18.59305556	Iceland
63.89777778	−18.05777778	Iceland
64.16555556	−17.74861111	Iceland
63.972315	−16.839655	Iceland
65.833333	−23.183333	Iceland
66.00507	−16.49788	Iceland
63.998806	−19.960389	Iceland
64.01325	−19.890556	Iceland
64.0695	−19.855722	Iceland
64.073472	−19.843444	Iceland
64.524833	−21.441028	Iceland
65.350694	−20.223972	Iceland
64.735431	−22.041482	Iceland
65.756178	−19.538124	Iceland
45.94	10.81	Italy
46.0137	11.0507	Italy
46.0346	11.018	Italy
46.0493	11.1957	Italy
46.027	10.9988	Italy
45.8709	10.7848	Italy
45.8843	10.8065	Italy
45.7526	10.8644	Italy
45.8266	10.9649	Italy
46.4741	11.7472	Italy
46.541592	11.71958	Italy
45.94001	10.81001	Italy
46.1279	11.2462	Italy
46.0166	11.2653	Italy
46.129173	11.360766	Italy
46.06118	11.2474	Italy
46.22	11.26	Italy
45.94002	10.81002	Italy
46.2398	11.236	Italy
46.2622	11.2709	Italy
46.275	11.2818	Italy
46.00564	11.084477	Italy
46.01371	11.05071	Italy
54.7976	25.3484	Lithuania
54.8921	24.1482	Lithuania
54.5729	24.6722	Lithuania
55.1526	25.8115	Lithuania
55.0938	26.0724	Lithuania
54.986	25.8061	Lithuania
54.6401	23.9624	Lithuania
69.9166	23.0001	NOR-ALT
70.1848	23.3726	NOR-ALT
59.381794	10.690076	NOR-SE
59.59135	5.1946	NOR-SW
59.032518	10.170094	NOR-SE
58.930892	9.595861	NOR-SE
69.5248	18.3752	NOR-TRO
59.262005	10.468178	NOR-SE
58.9323	9.53672	NOR-SE
58.783	9.2872	NOR-SE
68.7719	16.1878	NOR-LOF
59.0277	9.8772	NOR-SE
69.522	18.1593	NOR-TRO
58.494892	8.856994	NOR-SE
59.79135	5.522733	NOR-SW
59.772483	5.453716	NOR-SW
59.043257	9.689636	NOR-SE
63.341872	10.57585	NOR-CE
69.5359	20.3805	NOR-KOF
67.337324	15.564511	NOR-LOF
59.20753	10.93096	NOR-SE
59.168331	10.452577	NOR-SE
59.246213	10.405009	NOR-SE
63.33377	8.4582	NOR-CW
63.76167	11.60241	NOR-CE
63.80575	11.55992	NOR-CE
69.5274	20.3724	NOR-KOF
59.14931	11.125367	NOR-SE
63.902971	11.301787	NOR-CE
64.142987	11.276353	NOR-CE
62.913976	8.184851	NOR-CW
62.8904	8.5369	NOR-CW
69.47826	20.88754	NOR-KOF
58.91914	5.99499	NOR-SW
59.471639	6.25971	NOR-SW
62.28653	6.95902	NOR-CW
68.59102	14.867667	NOR-LOF
62.7346	7.0793	NOR-CW
59.325385	10.685341	NOR-SE
69.4451	20.949	NOR-KOF
59.382854	10.652894	NOR-SE
63.446121	11.349706	NOR-CE
63.692608	10.287903	NOR-CW
69.4256	20.9715	NOR-KOF
58.6469	9.1179	NOR-SE
63.132352	10.326228	NOR-CE
59.1406	6.0383	NOR-SW
68.80165	17.2212	NOR-LOF
69.90233	21.89166	NOR-ALT
70.0301	22.0653	NOR-ALT
69.9395	23.0964	NOR-ALT
70.0227	22.0173	NOR-ALT
70.1671	24.7561	NOR-POR
69.8537	25.051	NOR-POR
64.5	11.58	NOR-CE
69.694624	18.732784	NOR-TRO
69.691734	18.703207	NOR-TRO
69.530278	18.380833	NOR-TRO
69.9396	23.0971	NOR-ALT
69.783889	18.538889	NOR-TRO
70.0226	23.5595	NOR-ALT
70.0265	23.4942	NOR-ALT
70.0273	23.3874	NOR-ALT
70.0324	23.4012	NOR-ALT
70.1585	23.2781	NOR-ALT
38.7984	−9.38811	IBE-S
38.8	−9.383333333	IBE-S
46.77728792	23.55904282	Romania
62.268	33.982	Karelia
52.2869741	104.3050183	RUS-SIB
51.942186	85.9719355	RUS-SIB
NA	NA	RUS-SIB
58.370161	42.408667	Russia
58.391672	42.584989	Russia
58.514017	41.494	Russia
46.1904614	7.5449226	Switzerland
57.8954	−5.1613	Britain
51.484547	−1.026138	Britain
54.7331	−3.2105	Britain
50.2451	−3.77663	Britain
54.466609	−3.016579	Britain
54.6014	−3.3197	Britain
54.7261	−3.2176	Britain
51.8757	−4.9392	Britain
57.5667	−3.8833	Britain
54.565379	−3.136706	Britain
51.4845	−1.0261	Britain
", header = TRUE)

# Ensure Latitude and Longitude are numeric and filter NA
gps_data <- gps_data %>%
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>%
  filter(!is.na(Latitude) & !is.na(Longitude))

# 3. Define coordinates for all 15 FST regions
fst_coords <- data.frame(
  Region = fst_regions,
  Latitude = c(
    mean(gps_data$Latitude[gps_data$Region == "NOR-ALT"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "NOR-KOF"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "NOR-TRO"], na.rm = TRUE),
    64.5,  # Iceland
    mean(gps_data$Latitude[gps_data$Region == "NOR-CW"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "NOR-CE"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "NOR-SW"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "NOR-SE"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "Finland"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "Lithuania"], na.rm = TRUE),
    53.5,  # Britain
    mean(gps_data$Latitude[gps_data$Region == "Germany"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "Croatia"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "Italy"], na.rm = TRUE),
    mean(gps_data$Latitude[gps_data$Region == "IBE-NE"], na.rm = TRUE)
  ),
  Longitude = c(
    mean(gps_data$Longitude[gps_data$Region == "NOR-ALT"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "NOR-KOF"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "NOR-TRO"], na.rm = TRUE),
    -20.0,  # Iceland
    mean(gps_data$Longitude[gps_data$Region == "NOR-CW"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "NOR-CE"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "NOR-SW"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "NOR-SE"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "Finland"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "Lithuania"], na.rm = TRUE),
    -3.0,  # Britain
    mean(gps_data$Longitude[gps_data$Region == "Germany"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "Croatia"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "Italy"], na.rm = TRUE),
    mean(gps_data$Longitude[gps_data$Region == "IBE-NE"], na.rm = TRUE)
  )
)

# 4. Compute pairwise geographic distances
distance_matrix <- matrix(NA, n, n, dimnames = list(fst_regions, fst_regions))
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      distance_matrix[i, j] <- distHaversine(
        c(fst_coords$Longitude[i], fst_coords$Latitude[i]),
        c(fst_coords$Longitude[j], fst_coords$Latitude[j])
      ) / 1000  # Convert meters to kilometers
    }
  }
}

distances <- as.data.frame(as.table(distance_matrix)) %>%
  filter(!is.na(Freq)) %>%
  rename(Region1 = Var1, Region2 = Var2, Distance = Freq)

# 5. Convert FST matrix to long format
fst_long <- melt(as.matrix(fst_data), na.rm = TRUE) %>%
  rename(Region1 = Var1, Region2 = Var2, FST = value)

# 6. Merge FST and distance data
fst_dist_data <- full_join(fst_long, distances, by = c("Region1", "Region2")) %>%
  filter(!is.na(FST) & !is.na(Distance))

# 7. Define East/West classification for all 15 regions
region_direction <- data.frame(
  Region = fst_regions,
  Direction = c("east", "east", "west", "west", "west", "east", "west", "east", 
                "east", "east", "west", "west", "east", "west", "west")
)

# Add direction to fst_dist_data
fst_dist_data <- fst_dist_data %>%
  left_join(region_direction, by = c("Region1" = "Region")) %>%
  rename(Direction1 = Direction) %>%
  left_join(region_direction, by = c("Region2" = "Region")) %>%
  rename(Direction2 = Direction) %>%
  mutate(Pair_Type = case_when(
    Direction1 == "east" & Direction2 == "east" ~ "East-East",
    Direction1 == "west" & Direction2 == "west" ~ "West-West",
    TRUE ~ "East-West"
  ))

# 8. Compute correlations for each group
east_pairs <- fst_dist_data %>% filter(Pair_Type == "East-East")
west_pairs <- fst_dist_data %>% filter(Pair_Type == "West-West")
east_west_pairs <- fst_dist_data %>% filter(Pair_Type == "East-West")

cor_east <- cor.test(east_pairs$Distance, east_pairs$FST, method = "pearson")
cor_west <- cor.test(west_pairs$Distance, west_pairs$FST, method = "pearson")
cor_east_west <- cor.test(east_west_pairs$Distance, east_west_pairs$FST, method = "pearson")

# Extract correlation coefficients and p-values
cor_east_r <- round(cor_east$estimate, 3)
cor_east_p <- format.pval(cor_east$p.value, digits = 3)
cor_west_r <- round(cor_west$estimate, 3)
cor_west_p <- format.pval(cor_west$p.value, digits = 3)
cor_east_west_r <- round(cor_east_west$estimate, 3)
cor_east_west_p <- format.pval(cor_east_west$p.value, digits = 3)

# Create labels
east_label <- paste("East-East: r = ", cor_east_r, ", p = ", cor_east_p, sep = "")
west_label <- paste("West-West: r = ", cor_west_r, ", p = ", cor_west_p, sep = "")
east_west_label <- paste("East-West: r = ", cor_east_west_r, ", p = ", cor_east_west_p, sep = "")

# 9. Scatter plot with colored points and separate regression lines
scatter_plot <- ggplot(fst_dist_data, aes(x = Distance, y = FST, color = Pair_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", aes(color = Pair_Type), size = 1, se = FALSE) +
  scale_color_manual(values = c("East-East" = "red", "West-West" = "blue", "East-West" = "purple")) +
  scale_x_continuous(limits = c(0, 3500)) +  # X-axis from 0 to 3500 km
  scale_y_continuous(breaks = seq(0, 0.8, 0.1), limits = c(0, 0.8)) +  # Y-axis from 0 to 0.6, 0.1 increments
  labs(title = "Correlation between Geographic Distance and FST (All Regions)",
       x = "Geographic Distance (km)",
       y = expression(F[ST] * " (Genetic Distance)"),
       color = "Pair Type") +
  annotate("text", x = Inf, y = Inf, label = east_label, hjust = 1.1, vjust = 1.1, size = 4, color = "red") +
  annotate("text", x = Inf, y = Inf, label = west_label, hjust = 1.1, vjust = 2.5, size = 4, color = "blue") +
  annotate("text", x = Inf, y = Inf, label = east_west_label, hjust = 1.1, vjust = 4.0, size = 4, color = "purple") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the scatter plot
ggsave("fst_isolation_by_distance_all_regions_6_6.pdf", plot = scatter_plot, width = 6, height = 6, dpi = 300)

# 10. Statistical test for differences in FST between groups
kruskal_test <- kruskal.test(FST ~ Pair_Type, data = fst_dist_data)
pairwise_wilcox <- pairwise.wilcox.test(fst_dist_data$FST, fst_dist_data$Pair_Type, 
                                        p.adjust.method = "bonferroni")

# Print test results
cat("Kruskal-Wallis test for overall differences:\n")
print(kruskal_test)
cat("Pairwise Wilcoxon tests with Bonferroni correction:\n")
print(pairwise_wilcox)

# Prepare comparisons for annotation
comparisons <- list(c("East-East", "West-West"), 
                    c("East-West", "West-West"), 
                    c("East-East", "East-West"))

# 11. Boxplot of FST values by group with significance
box_plot <- ggplot(fst_dist_data, aes(x = Pair_Type, y = FST, fill = Pair_Type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("East-East" = "red", "West-West" = "blue", "East-West" = "purple")) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.1), limits = c(0, 0.8)) +  # Y-axis from 0 to 0.6, 0.1 increments
  labs(title = "FST Distribution by Pair Type (All Regions)",
       x = "Pair Type",
       y = expression(F[ST] * " (Genetic Distance)"),
       fill = "Pair Type") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 0.58, label.x = 1.5) +  # Overall Kruskal p-value
  stat_compare_means(aes(group = Pair_Type), method = "wilcox.test", 
                     label = "p.signif", label.y = c(0.5, 0.52, 0.54), 
                     comparisons = comparisons)  # Pairwise comparisons

# Save the boxplot
ggsave("fst_boxplot_all_regions_with_significance.pdf", plot = box_plot, width = 6, height = 6, dpi = 300)

# 12. Print correlations for reference
cat("Eastern pairs correlation:\n")
print(cor_east)
cat("Western pairs correlation:\n")
print(cor_west)
cat("East-West pairs correlation:\n")
print(cor_east_west)

# 13. Verify the number of pairs
cat("Total unique regions:", length(unique(c(fst_dist_data$Region1, fst_dist_data$Region2))), "\n")
cat("Total pairs in dataset:", nrow(fst_dist_data), "\n")
cat("Expected pairs for 15 regions:", choose(15, 2), "\n")
cat("Number of East-East pairs:", nrow(east_pairs), "\n")
cat("Number of West-West pairs:", nrow(west_pairs), "\n")
cat("Number of East-West pairs:", nrow(east_west_pairs), "\n")